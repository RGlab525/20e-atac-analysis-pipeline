# =============================================================================
# Supplementary sample QC (12 libraries; Editor 1, Reviewer 3 Minor 3, Reviewer 1 Minor 1)
#   Input BAMs are the filtered analysis BAMs (unmapped, secondary/supplementary, QC-fail and MAPQ < 20 reads
#   already removed; duplicates were not marked or removed), so raw mapping rate and duplicate % cannot be
#   computed here.
#   - filtered fragments (proper pairs), chrM fragments and %
#   - FRiP: fraction of nuclear fragments overlapping the 14,003 consensus peaks
#   - TSS enrichment (ENCODE style): Tn5 insertion sites (+4/-5) aggregated over TSS +/- 2 kb (org.Dm.eg.db
#     CHRLOC, strand-aware), normalized to the mean of the outer 100 bp on each side, maximum of the 11-bp
#     running mean within +/- 50 bp of the TSS
#   - fragment length: % < 100 bp, % 180-247 bp (mono-nucleosome), median
# Input : input_bam/input_atac_rep1_{t}.sorted.bam (rep1), input_bam/input_atac_rep2_{t}.sorted.bam (rep2), input_consensus_peaks.bed
# Output: SuppFig/output_ATAC_sample_QC.csv, SuppFig/output_TSS_profiles.rds
# Source: 261016_revision/08_QC_peak_validity_fragment/scripts/atac_sample_qc.R (J.Y.H.)
# Changes from the original code: paths via config/config.R; output names without the date stamp.
# =============================================================================

suppressMessages({
    library(Rsamtools)
    library(GenomicRanges)
    library(org.Dm.eg.db)
})

source("config/config.R")
proj   <- PROJECT_DIR
outdir <- normalizePath(out_dir("SuppFig"))
proj   <- normalizePath(proj)
setwd(file.path(proj, "input_bam"))

# ----- peaks -----
pk <- read.table(file.path(proj, "input_consensus_peaks.bed"), sep = "\t")
peaks <- GRanges(sub("^chr", "", pk$V1), IRanges(pk$V2 + 1, pk$V3))

# ----- TSS (org.Dm.eg.db CHRLOC: negative = minus strand) -----
loc  <- toTable(org.Dm.egCHRLOC)
locE <- toTable(org.Dm.egCHRLOCEND)
names(loc)[2] <- "start"; names(locE)[2] <- "end"
tss <- merge(loc, locE, by = c("gene_id", "Chromosome"))
tss <- tss[sign(tss$start) == sign(tss$end), ]
tss$strand <- ifelse(tss$start > 0, "+", "-")
tss$pos    <- ifelse(tss$strand == "+", abs(tss$start), abs(tss$end))
tss <- unique(tss[tss$Chromosome %in% c("2L", "2R", "3L", "3R", "4", "X"), c("Chromosome", "pos", "strand")])
cat("TSS:", nrow(tss), "\n")

W <- 2000
tss.enrichment <- function(chr, ins) {
    prof <- numeric(2 * W + 1)
    for (ch in unique(tss$Chromosome)) {
        t <- tss[tss$Chromosome == ch, ]
        x <- sort(ins[chr == ch])
        if (!length(x)) next
        for (s in c("+", "-")) {
            tp <- t$pos[t$strand == s]
            for (p in tp) {
                lo <- findInterval(p - W - 1, x) + 1; hi <- findInterval(p + W, x)
                if (hi < lo) next
                d <- x[lo:hi] - p
                if (s == "-") d <- -d
                tb <- tabulate(d + W + 1, nbins = 2 * W + 1)
                prof <- prof + tb
            }
        }
    }
    flank <- mean(c(prof[1:100], prof[(2 * W - 98):(2 * W + 1)]))
    norm  <- prof / flank
    sm    <- stats::filter(norm, rep(1 / 11, 11), sides = 2)
    list(score = max(sm[(W + 1 - 50):(W + 1 + 50)], na.rm = TRUE), profile = norm)
}

qc.one <- function(bam, sample, rep) {
    message("== ", sample, " ", rep, " : ", bam)
    p <- ScanBamParam(what = c("rname", "pos", "isize"),
                      flag = scanBamFlag(isProperPair = TRUE, isFirstMateRead = TRUE, isUnmappedQuery = FALSE))
    x <- scanBam(bam, param = p)[[1]]
    keep <- !is.na(x$isize) & x$isize != 0
    chr <- as.character(x$rname[keep]); pos <- x$pos[keep]; isz <- x$isize[keep]
    fs  <- ifelse(isz > 0, pos, pos + isz); len <- abs(isz)
    fe  <- fs + len - 1
    n.all  <- length(fs)
    is.m   <- chr %in% c("M", "chrM", "MT")
    n.m    <- sum(is.m)
    chr <- chr[!is.m]; fs <- fs[!is.m]; fe <- fe[!is.m]; len <- len[!is.m]
    fr  <- GRanges(chr, IRanges(fs, fe))
    in.peak <- overlapsAny(fr, peaks)
    ins.chr <- c(chr, chr); ins <- c(fs + 4, fe - 5)
    te <- tss.enrichment(ins.chr, ins)
    res <- data.frame(sample = sample, replicate = rep,
               filtered_fragments = n.all,
               chrM_fragments = n.m, chrM_pct = round(100 * n.m / n.all, 2),
               nuclear_fragments = length(fs),
               FRiP_pct = round(100 * mean(in.peak), 2),
               TSS_enrichment = round(te$score, 2),
               frag_lt100_pct = round(100 * mean(len < 100), 2),
               frag_mono_180_247_pct = round(100 * mean(len >= 180 & len <= 247), 2),
               median_fragment_length = median(len))
    attr(res, "profile") <- te$profile
    res
}

tp <- c("0h", "1h", "3h", "6h", "12h", "24h")
jobs <- rbind(data.frame(t = tp, rep = "rep1", bam = sprintf("input_atac_rep1_%s.sorted.bam", tp)),
              data.frame(t = tp, rep = "rep2", bam = sprintf("input_atac_rep2_%s.sorted.bam", tp)))
out <- list(); profs <- list()
for (i in seq_len(nrow(jobs))) {
    r <- qc.one(jobs$bam[i], jobs$t[i], jobs$rep[i])
    profs[[paste(jobs$t[i], jobs$rep[i])]] <- attr(r, "profile")
    attr(r, "profile") <- NULL
    out[[i]] <- r
    print(r, row.names = FALSE)
    write.csv(do.call(rbind, out), file.path(outdir, "output_ATAC_sample_QC.csv"), row.names = FALSE)
}
qc <- do.call(rbind, out)
qc$sample <- factor(qc$sample, levels = tp); qc <- qc[order(qc$sample, qc$replicate), ]
write.csv(qc, file.path(outdir, "output_ATAC_sample_QC.csv"), row.names = FALSE)
saveRDS(profs, file.path(outdir, "output_TSS_profiles.rds"))
cat("\n== Summary ==\n"); print(qc, row.names = FALSE)
