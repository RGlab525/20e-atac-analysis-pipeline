# =============================================================================
# Fig 2G (new box-plot panel), step 1: windows in which the mean ATAC signal is measured
#   TSS peaks      : -425/+165 bp around the representative TSS of each of the 4,774 TSS peaks (strand-aware;
#                    the TSS-proximal window defined in Supplementary Fig. S3A)
#   distal peaks   : peak centre +/-200 bp for the 9,229 distal peaks
#   Intergenic values are not made here: they come from the Supplementary Fig. S2A matrix
#   (SuppFig_S2_02_intergenic_genebody_fingerprint.sh -> SuppFig_S2_03_intergenic_genebody_values.R).
#
# Input:  $ATAC_OUT/Fig2AB_regions/output_A_TSS_peaks_4774.bed, output_B_distal_peaks_9229.bed   (Fig2AB_01_region_beds.R)
#         input_bam/input_atac_rep1_0h.sorted.bam   (only its header, for chromosome lengths)
# Output: $ATAC_OUT/Fig2G_windows/windows_TSS_peaks.bed, windows_ENH_distal9229.bed
#
# Source: 261016_revision/03_Fig2AB_heatmap_boxplot/scripts/01_make_boxplot_windows.R
#   (the windows of the original Fig 2A/B sets, made by the same script for comparison, are not repeated here)
# =============================================================================
source("config/config.R")
reg <- file.path(OUT_ROOT, "Fig2AB_regions")
out <- out_dir("Fig2G_windows")
chrs <- paste0("chr", c("2L", "2R", "3L", "3R", "4", "X", "Y", "M"))   # all chromosomes, so the 9,229 distal peaks match Figs 3-6

tss.win <- function(d) {                       # d: BED6 with 1-bp TSS records
    d <- d[d$V1 %in% chrs, ]; pos <- d$V3      # 1-based TSS position
    s <- ifelse(d$V6 == "-", pos - 165, pos - 425); e <- ifelse(d$V6 == "-", pos + 425, pos + 165)
    unique(data.frame(chr = d$V1, start = s - 1, end = e))
}
ctr.win <- function(d, half = 200) {
    d <- d[d$V1 %in% chrs, ]; c <- (d$V2 + d$V3) %/% 2
    unique(data.frame(chr = d$V1, start = c - half, end = c + half))
}
rd <- function(f) read.table(f, sep = "\t", quote = "", comment.char = "", fill = TRUE)
w <- list(
    TSS_peaks      = tss.win(rd(file.path(reg, "output_A_TSS_peaks_4774.bed"))),
    ENH_distal9229 = ctr.win(rd(file.path(reg, "output_B_distal_peaks_9229.bed"))))
len <- Rsamtools::scanBamHeader(inp("input_bam/input_atac_rep1_0h.sorted.bam"))[[1]]$targets   # chromosome lengths (same dm6 as the bigWigs)
names(len) <- paste0("chr", names(len))
for (n in names(w)) {
    x <- w[[n]]; x$start <- pmax(x$start, 0); x$end <- pmin(x$end, len[x$chr])          # clip windows at chromosome ends (chrM)
    x <- unique(x[x$end > x$start, ])
    x <- x[order(x$chr, x$start), ]; x$name <- paste0(n, "_", seq_len(nrow(x)))
    write.table(x, file.path(out, paste0("output_windows_", n, ".bed")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat(sprintf("%-18s %6d windows\n", n, nrow(x)))
}
