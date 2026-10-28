# =============================================================================
# R3 Major 8 (E1/E2) and R3 Major 7 / R2-5: overlap of our ATAC peaks with public S2 STARR-seq enhancers
#   (Shlyueva et al. 2014, GSE47691; untreated and 20E 41 uM 24 h, 2 replicates each)
#   STARR-seq peak file = 1-bp summit (dm3) with enrichment and binomial P; summit +-250 bp is used as the enhancer interval
#   main cutoff P <= 0.001 (also reported: P <= 0.001 and enrichment >= 3); DHS-seq peaks lifted the same way (overlap only)
#   dm3 -> dm6 with UCSC dm3ToDm6.over.chain.gz (rtracklayer::liftOver)
#   genome-wide: distal peaks (9,229) increased / decreased / unchanged at 3 h vs STARR-seq; expectation from 100
#   within-chromosome shuffles of the distal peaks (set.seed(1)); E1 rank among STARR-seq summits per sample
# Input:  public_data/GSE47691_STARR_DHS/ (STARRseq_GSE47691_00_download.sh),
#         input_peak_annotations.xlsx (sheet 2), input_consensus_peaks.bed
# Output: $ATAC_OUT/STARRseq_GSE47691/: output_CR43626_locus_STARRseq_overlap.csv, output_CR43626_locus_DHSseq_overlap.csv,
#         output_CR43626_region_all_STARR_summits_dm6.csv (track input of Supplementary Fig S7), output_genomewide_peaks_STARRseq_overlap.csv,
#         output_genomewide_peaks_DHSseq_overlap.csv, output_distal_DA3h_vs_STARR_fisher.csv,
#         output_E1_STARR_rank_per_sample.csv, output_lifted_dm6.rds
# Source: 261016_revision/15_R3_STARRseq_E1E2/scripts/01_starr_overlap.R (paths only changed)
# =============================================================================
source("config/config.R")
suppressPackageStartupMessages({library(rtracklayer); library(GenomicRanges); library(readxl)})
dat <- inp("public_data/GSE47691_STARR_DHS"); res <- out_dir("STARRseq_GSE47691")
P_CUT <- 1e-3; ENR_CUT <- 3; HALF <- 250

chain <- import.chain(gzcon_tmp <- {f <- tempfile(fileext = ".chain"); system2("gunzip", c("-c", shQuote(file.path(dat, "input_dm3ToDm6.over.chain.gz"))), stdout = f); f})

# ---- STARR-seq ----
starr_files <- c(noEcd_rep1 = "input_GSM1357045_S2_STARRseq_no_ecd_rep1", noEcd_rep2 = "input_GSM1357046_S2_STARRseq_no_ecd_rep2",
                 ecd_rep1   = "input_GSM1357047_S2_STARRseq_with_ecd_rep1", ecd_rep2 = "input_GSM1357048_S2_STARRseq_with_ecd_rep2")
lift_summits <- function(f) {
  x <- read.table(file.path(dat, paste0(f, ".peaks.txt.gz")), col.names = c("chr", "summit", "enrichment", "P"))
  g <- GRanges(x$chr, IRanges(x$summit, width = 1), enrichment = x$enrichment, P = x$P, summit_dm3 = x$summit)
  l <- liftOver(g, chain); ok <- lengths(l) == 1
  out <- unlist(l[ok]); out$summit_dm6 <- start(out)
  attr(out, "n") <- c(total = length(g), lifted = sum(ok)); out
}
starr <- lapply(starr_files, lift_summits)
win   <- lapply(starr, function(s) resize(s, 2 * HALF + 1, fix = "center"))
cat("STARR-seq summits (total / lifted to dm6 / P<=0.001 / P<=0.001 & enr>=3):\n")
for (k in names(starr)) cat(sprintf("  %-10s %5d / %5d / %5d / %5d\n", k, attr(starr[[k]], "n")[1], attr(starr[[k]], "n")[2],
                                    sum(starr[[k]]$P <= P_CUT), sum(starr[[k]]$P <= P_CUT & starr[[k]]$enrichment >= ENR_CUT)))

# ---- DHS-seq (MACS intervals) ----
dhs_files <- c(noEcd_rep1 = "input_GSM1154837_S2_DHSseq_no_ecd_rep1", noEcd_rep2 = "input_GSM1154838_S2_DHSseq_no_ecd_rep2",
               ecd_rep1   = "input_GSM1154839_S2_DHSseq_with_ecd_rep1", ecd_rep2 = "input_GSM1154840_S2_DHSseq_with_ecd_rep2")
dhs <- lapply(dhs_files, function(f) {
  x <- read.table(file.path(dat, paste0(f, ".peaks.txt.gz")), header = TRUE, comment.char = "#", check.names = FALSE)
  g <- GRanges(x$chr, IRanges(x$start, x$end), fold = x$fold_enrichment)
  l <- liftOver(g, chain); keep <- lengths(l) > 0
  r <- range(l[keep]); ok <- lengths(r) == 1      # keep intervals lifted as one piece
  out <- unlist(r[ok]); out$fold <- g$fold[keep][ok]; out
})

# ---- our peaks ----
ann <- as.data.frame(read_excel(inp("input_peak_annotations.xlsx"), sheet = 2))
pk  <- GRanges(ann$chr, IRanges(ann$start, ann$end), peak_id = ann$peak_id, TSS = ann$TSS,
               lfc3 = ann$log2FoldChange_3h, padj3 = ann$padj_3h)
# check that the table matches the BED file (the table keeps the BED numbers)
bed <- read.table(inp("input_consensus_peaks.bed"))
stopifnot(all(ann$peak_id %in% bed$V4), nrow(bed) == nrow(ann))

# ---- locus: E1, P, E2 + other CR43626-linked distal peaks ----
loc <- data.frame(region = c("peak6199", "E1", "P", "E2", "peak6203", "peak6204"),
                  peak_id = c("peak6199", "peak6200", "peak6201", "peak6202", "peak6203", "peak6204"))
lg <- pk[match(loc$peak_id, pk$peak_id)]
loc <- cbind(loc, coord = sprintf("chr3L:%s-%s", format(start(lg), big.mark = ","), format(end(lg), big.mark = ",")),
             ATAC_3h_log2FC = round(lg$lfc3, 2), ATAC_3h_padj = signif(lg$padj3, 2))
rows <- list()
for (i in seq_along(lg)) for (k in names(starr)) {
  s <- starr[[k]]; d <- distance(lg[i], s); d[is.na(d)] <- Inf
  h <- which(overlapsAny(win[[k]], lg[i]))                     # summit +-250 bp overlaps the peak
  hs <- h[s$P[h] <= P_CUT]
  j <- if (length(hs)) hs[which.max(s$enrichment[hs])] else if (length(h)) h[which.max(s$enrichment[h])] else NA
  nn <- which.min(ifelse(s$P <= P_CUT, d, Inf))                 # nearest P <= 0.001 summit
  rows[[length(rows) + 1]] <- data.frame(region = loc$region[i], sample = k,
    overlap_P001 = length(hs) > 0, n_summits_overlap_any = length(h),
    best_summit_dm6 = if (is.na(j)) NA else s$summit_dm6[j], summit_in_peak = if (is.na(j)) NA else d[j] == 0,
    enrichment = if (is.na(j)) NA else round(s$enrichment[j], 1), P = if (is.na(j)) NA else signif(s$P[j], 2),
    nearest_P001_summit_dm6 = s$summit_dm6[nn], nearest_P001_dist_bp = d[nn], nearest_P001_enrichment = round(s$enrichment[nn], 1))
}
locus <- merge(loc, do.call(rbind, rows), by = "region", sort = FALSE)
locus <- locus[order(match(locus$region, loc$region), match(locus$sample, names(starr))), ]
for (k in names(dhs)) {
  o <- findOverlaps(lg, dhs[[k]])
  f <- tapply(dhs[[k]]$fold[subjectHits(o)], queryHits(o), max)
  loc[[paste0("DHS_", k)]] <- ifelse(seq_along(lg) %in% queryHits(o), paste0("yes (FE ", f[as.character(seq_along(lg))], ")"), "no")
}
write.csv(locus, file.path(res, "output_CR43626_locus_STARRseq_overlap.csv"), row.names = FALSE)
write.csv(loc,   file.path(res, "output_CR43626_locus_DHSseq_overlap.csv"),   row.names = FALSE)
# all nearby summits (for the figure / checking)
near <- GRanges("chr3L", IRanges(3235000, 3262000))
nb <- do.call(rbind, lapply(names(starr), function(k) { s <- starr[[k]][overlapsAny(starr[[k]], near)]
  if (!length(s)) return(NULL); data.frame(sample = k, summit_dm3 = s$summit_dm3, summit_dm6 = s$summit_dm6, enrichment = round(s$enrichment, 2), P = signif(s$P, 3)) }))
write.csv(nb[order(nb$sample, nb$summit_dm6), ], file.path(res, "output_CR43626_region_all_STARR_summits_dm6.csv"), row.names = FALSE)

# ---- genome-wide ----
sets <- list(distal_all = pk[pk$TSS == "X"], TSS_all = pk[pk$TSS == "O"])
d <- sets$distal_all
sets$distal_up3h   <- d[!is.na(d$padj3) & d$padj3 < 0.05 & d$lfc3 > 0]
sets$distal_down3h <- d[!is.na(d$padj3) & d$padj3 < 0.05 & d$lfc3 < 0]
sets$distal_ns3h   <- d[is.na(d$padj3) | d$padj3 >= 0.05]
stopifnot(length(sets$distal_all) == 9229)
sel <- function(k, strict = FALSE) { s <- win[[k]]; s[s$P <= P_CUT & (!strict | s$enrichment >= ENR_CUT)] }
cond <- list(noEcd_any = c("noEcd_rep1", "noEcd_rep2"), ecd_any = c("ecd_rep1", "ecd_rep2"),
             any_of_4 = names(starr))
hit <- function(g, ks, strict = FALSE) Reduce(`|`, lapply(ks, function(k) overlapsAny(g, sel(k, strict))))
both <- function(g, ks, strict = FALSE) Reduce(`&`, lapply(ks, function(k) overlapsAny(g, sel(k, strict))))

# expectation: distal peaks shuffled within chromosomes (100 times)
set.seed(1)
seqlen <- tapply(end(pk), as.character(seqnames(pk)), max) + 5000
shuffle <- function(g) { L <- seqlen[as.character(seqnames(g))]
  st <- floor(runif(length(g), 1, L - width(g))); GRanges(seqnames(g), IRanges(st, width = width(g))) }
perm <- replicate(100, { r <- shuffle(sets$distal_all); c(any_of_4 = mean(hit(r, cond$any_of_4)),
  noEcd_any = mean(hit(r, cond$noEcd_any)), ecd_any = mean(hit(r, cond$ecd_any))) })

gw <- do.call(rbind, lapply(names(sets), function(sn) { g <- sets[[sn]]
  do.call(rbind, lapply(names(cond), function(cn) data.frame(set = sn, n_peaks = length(g), STARR = cn,
    n_overlap_P001_anyRep = sum(hit(g, cond[[cn]])), pct_P001_anyRep = round(100 * mean(hit(g, cond[[cn]])), 1),
    n_overlap_P001_bothReps = if (cn == "any_of_4") NA else sum(both(g, cond[[cn]])),
    n_overlap_P001_enr3_anyRep = sum(hit(g, cond[[cn]], TRUE)), pct_P001_enr3_anyRep = round(100 * mean(hit(g, cond[[cn]], TRUE)), 1),
    pct_expected_shuffled_distal = if (sn == "distal_all") round(100 * mean(perm[cn, ]), 1) else NA)))
}))
write.csv(gw, file.path(res, "output_genomewide_peaks_STARRseq_overlap.csv"), row.names = FALSE)

# DHS genome-wide (reference)
gd <- do.call(rbind, lapply(c("distal_all", "TSS_all"), function(sn) data.frame(set = sn,
  pct_DHS_any_of_4 = round(100 * mean(Reduce(`|`, lapply(dhs, function(x) overlapsAny(sets[[sn]], x)))), 1))))
write.csv(gd, file.path(res, "output_genomewide_peaks_DHSseq_overlap.csv"), row.names = FALSE)

# Fisher test of 3 h DA vs STARR overlap (within distal peaks)
f_tab <- function(a, b) { x <- matrix(c(sum(hit(a, cond$any_of_4)), length(a) - sum(hit(a, cond$any_of_4)),
                                       sum(hit(b, cond$any_of_4)), length(b) - sum(hit(b, cond$any_of_4))), 2)
  ft <- fisher.test(x); c(OR = unname(round(ft$estimate, 2)), P = signif(ft$p.value, 2)) }
fis <- rbind(up_vs_ns = f_tab(sets$distal_up3h, sets$distal_ns3h), down_vs_ns = f_tab(sets$distal_down3h, sets$distal_ns3h))
f_tab2 <- function(a, b, ks) { x <- matrix(c(sum(hit(a, ks)), length(a) - sum(hit(a, ks)), sum(hit(b, ks)), length(b) - sum(hit(b, ks))), 2)
  ft <- fisher.test(x); c(OR = unname(round(ft$estimate, 2)), P = signif(ft$p.value, 2)) }
fis <- rbind(fis, up_vs_ns_ecdSTARR = f_tab2(sets$distal_up3h, sets$distal_ns3h, cond$ecd_any),
             down_vs_ns_noEcdSTARR = f_tab2(sets$distal_down3h, sets$distal_ns3h, cond$noEcd_any))
fis <- cbind(STARR = c("any_of_4", "any_of_4", "ecd_any", "noEcd_any"), fis)
write.csv(fis, file.path(res, "output_distal_DA3h_vs_STARR_fisher.csv"))

# rank of the E1 summit per sample (among P <= 0.001 summits; ranks because enrichment scales differ between samples)
e1 <- lg[2]
rk <- do.call(rbind, lapply(names(starr), function(k) { s <- starr[[k]][starr[[k]]$P <= P_CUT]
  h <- which(overlapsAny(resize(s, 2 * HALF + 1, fix = "center"), e1)); if (!length(h)) return(NULL)
  e <- max(s$enrichment[h]); r <- sum(s$enrichment > e) + 1
  data.frame(sample = k, E1_enrichment = round(e, 1), rank = r, n_P001 = length(s), top_pct = round(100 * r / length(s), 1),
             median_enrichment_P001 = round(median(s$enrichment), 1)) }))
write.csv(rk, file.path(res, "output_E1_STARR_rank_per_sample.csv"), row.names = FALSE)

options(width = 200)
cat("\n== E1 rank ==\n"); print(rk, row.names = FALSE)
cat("\n== locus (STARR) ==\n"); print(locus[, c("region", "sample", "overlap_P001", "n_summits_overlap_any", "best_summit_dm6", "enrichment", "P", "nearest_P001_dist_bp", "nearest_P001_enrichment")], row.names = FALSE)
cat("\n== locus (DHS) ==\n"); print(loc, row.names = FALSE)
cat("\n== genome-wide ==\n"); print(gw, row.names = FALSE)
cat("\n"); print(gd, row.names = FALSE); cat("\n"); print(fis)
cat("\nshuffled expectation (mean %):\n"); print(round(100 * rowMeans(perm), 1))
saveRDS(list(starr = starr, dhs = dhs), file.path(res, "output_lifted_dm6.rds"))
