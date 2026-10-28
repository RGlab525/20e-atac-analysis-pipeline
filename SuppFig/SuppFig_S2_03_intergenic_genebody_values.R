# =============================================================================
# Supplementary Fig. S2, step 3: mean signal per intergenic region and per gene body (input of the box plots)
#   For each region and time point: mean of the 300 bins (3,000 bp / 10 bp) of the computeMatrix matrix.
#   Rows with the same coordinates (several transcripts) are one locus and are kept once.
#   "noPeak" sets: regions that do not overlap any consensus peak (background-only check).
#
# Input:  $ATAC_OUT/SuppFig_S2/output_Intergenic_S2A_down3000.mat.gz, output_GeneBody_TSSplus1kb_TES.mat.gz (step 2)
#         input_consensus_peaks.bed
# Output ($ATAC_OUT/SuppFig_S2/boxplot/, same format as multiBigwigSummary --outRawCounts):
#   Intergenic_mean_signal.tab, Intergenic_noPeak_mean_signal.tab, GeneBody_mean_signal.tab, GeneBody_noPeak_mean_signal.tab
#
# Source: 261016_revision/03_Fig2AB_heatmap_boxplot/scripts/01b_intergenic_from_original_S2A_matrix.R and
#   01d_genebody_from_matrix.R (merged; same calculation)
# =============================================================================
source("config/config.R")
suppressMessages(library(GenomicRanges)); options(warn = -1)
M <- file.path(OUT_ROOT, "SuppFig_S2")
O <- out_dir("SuppFig_S2/boxplot")
pk <- read.table(inp("input_consensus_peaks.bed"), sep = "\t")
pkg <- GRanges(pk$V1, IRanges(pk$V2 + 1, pk$V3))

region.means <- function(matfile, name) {
    m <- read.table(file.path(M, matfile), sep = "\t", skip = 1, quote = "", comment.char = "", na.strings = "nan")
    stopifnot(ncol(m) == 6 + 6 * 300)                                   # 6 time points x 300 bins
    vals <- sapply(seq_along(TIMES), function(i) rowMeans(as.matrix(m[, 6 + ((i - 1) * 300 + 1):(i * 300)]), na.rm = TRUE))
    colnames(vals) <- TIMES
    d <- data.frame(chr = paste0("chr", sub("^chr", "", m$V1)), start = m$V2, end = m$V3, vals, check.names = FALSE)
    cat(name, "- matrix rows:", nrow(d))
    d <- d[!duplicated(d[, c("chr", "start", "end")]), ]
    hasPeak <- overlapsAny(GRanges(d$chr, IRanges(d$start + 1, d$end)), pkg)
    cat("; unique loci:", nrow(d), "; overlapping a consensus peak:", sum(hasPeak), "; without peak:", sum(!hasPeak), "\n")
    wr <- function(x, f) { names(x)[1:3] <- c("#'chr'", "'start'", "'end'"); names(x)[4:9] <- paste0("'", TIMES, "'")
        write.table(x, file.path(O, f), sep = "\t", quote = FALSE, row.names = FALSE) }
    wr(d, paste0("output_", name, "_mean_signal.tab")); wr(d[!hasPeak, ], paste0("output_", name, "_noPeak_mean_signal.tab"))
    cat(sprintf("  median 0h %.3f -> 24h %.3f\n", median(d$`0h`), median(d$`24h`)))
}
region.means("output_Intergenic_S2A_down3000.mat.gz", "Intergenic")
region.means("output_GeneBody_TSSplus1kb_TES.mat.gz", "GeneBody")
