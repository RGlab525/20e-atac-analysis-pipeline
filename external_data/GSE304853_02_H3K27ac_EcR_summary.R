# =============================================================================
# R3 Major 8: H3K27ac / EcR change (GSE304853, 20E 0.3 uM 1 h vs DMSO) at E1, P, E2 of the CR43626 / miR-282 locus,
#   their percentile among all distal peaks, the Fig 1C/D boxed peaks as positive controls, and the genome-wide relation
#   between ATAC DAR groups and H3K27ac change. Numbers are reported in the response letter (R3 Major 8).
#   change = log2((20E + 0.1) / (DMSO + 0.1)), per replicate and for the replicate mean
#
# Input:  $ATAC_OUT/GSE304853/output_peaks_H3K27ac_EcR_20E_vs_DMSO.csv (GSE304853_01_H3K27ac_EcR_at_peaks.py)
#         input_peak_annotations.xlsx (sheet 2: ATAC DESeq2 1 h and 3 h per peak)
#         BOX_PEAKS: per-peak statistics of the Fig 1C/D boxed peaks (peak, gene, ...; made with the Fig 1C/D
#                    track figure, whose drawing code is not included), default input_box_peaks.csv
# Output: $ATAC_OUT/GSE304853/output_peaks_H3K27ac_EcR_with_ATAC.csv; summary tables printed to the console
# Source: 261016_revision/12_R3_Major8_H3K27ac/scripts/02_summary.R (paths only changed)
# =============================================================================
source("config/config.R")
suppressMessages(library(readxl)); options(warn = -1)
BOX_PEAKS <- Sys.getenv("BOX_PEAKS", "input_box_peaks.csv")
B <- out_dir("GSE304853")
x <- read.csv(file.path(B, "output_peaks_H3K27ac_EcR_20E_vs_DMSO.csv"))
a <- as.data.frame(read_excel(inp("input_peak_annotations.xlsx"), sheet = 2))
x <- cbind(x, a[match(x$peak, a$peak_id), c("TSS", "log2FoldChange_1h", "padj_1h", "log2FoldChange_3h", "padj_3h")])
lf <- function(t, c) log2((t + 0.1) / (c + 0.1))
x$K27_lfc   <- lf((x$K27_20E_rep1 + x$K27_20E_rep2) / 2, (x$K27_DMSO_rep1 + x$K27_DMSO_rep2) / 2)
x$K27_lfc_r1 <- lf(x$K27_20E_rep1, x$K27_DMSO_rep1); x$K27_lfc_r2 <- lf(x$K27_20E_rep2, x$K27_DMSO_rep2)
x$EcR_lfc   <- lf((x$EcR_20E_rep1 + x$EcR_20E_rep2) / 2, (x$EcR_DMSO_rep1 + x$EcR_DMSO_rep2) / 2)
x$EcR_lfc_r1 <- lf(x$EcR_20E_rep1, x$EcR_DMSO_rep1); x$EcR_lfc_r2 <- lf(x$EcR_20E_rep2, x$EcR_DMSO_rep2)
x$K27_DMSO <- (x$K27_DMSO_rep1 + x$K27_DMSO_rep2) / 2; x$K27_20E <- (x$K27_20E_rep1 + x$K27_20E_rep2) / 2
x$EcR_DMSO <- (x$EcR_DMSO_rep1 + x$EcR_DMSO_rep2) / 2; x$EcR_20E <- (x$EcR_20E_rep1 + x$EcR_20E_rep2) / 2
x$class <- ifelse(x$TSS == "O", "TSS", "distal")
fc <- log2(1.5)
x$atac1 <- ifelse(x$padj_1h < 0.05 & x$log2FoldChange_1h > fc, "up", ifelse(x$padj_1h < 0.05 & x$log2FoldChange_1h < -fc, "down", "unchanged"))
x$atac3 <- ifelse(x$padj_3h < 0.05 & x$log2FoldChange_3h > fc, "up", ifelse(x$padj_3h < 0.05 & x$log2FoldChange_3h < -fc, "down", "unchanged"))
x$atac1[is.na(x$atac1)] <- "unchanged"; x$atac3[is.na(x$atac3)] <- "unchanged"
write.csv(x, file.path(B, "output_peaks_H3K27ac_EcR_with_ATAC.csv"), row.names = FALSE)

# 1) E1, P, E2
key <- data.frame(region = c("E1", "P", "E2"), peak = c("peak6200", "peak6201", "peak6202"))
d <- x[x$class == "distal", ]
pct <- function(v, ref) round(100 * mean(ref <= v), 1)
k <- merge(key, x, by = "peak")[match(key$peak, merge(key, x, by = "peak")$peak), ]
cat("== E1 / P / E2 (ATAC = submitted DESeq2; H3K27ac, EcR = GSE304853, 20E 0.3 uM 1 h vs DMSO)\n")
print(data.frame(region = k$region, peak = k$peak, coord = paste0(k$chr, ":", k$start, "-", k$end),
                 ATAC_1h = round(k$log2FoldChange_1h, 2), ATAC_3h = round(k$log2FoldChange_3h, 2), padj_3h = signif(k$padj_3h, 2),
                 K27_DMSO = round(k$K27_DMSO, 2), K27_20E = round(k$K27_20E, 2), K27_lfc = round(k$K27_lfc, 2),
                 K27_rep1 = round(k$K27_lfc_r1, 2), K27_rep2 = round(k$K27_lfc_r2, 2),
                 K27_pctile_distal = sapply(k$K27_lfc, pct, d$K27_lfc),
                 EcR_DMSO = round(k$EcR_DMSO, 2), EcR_20E = round(k$EcR_20E, 2), EcR_lfc = round(k$EcR_lfc, 2),
                 EcR_rep1 = round(k$EcR_lfc_r1, 2), EcR_rep2 = round(k$EcR_lfc_r2, 2)), row.names = FALSE)
cat("distal peaks: H3K27ac DMSO signal median", round(median(d$K27_DMSO), 2), " 90th pct", round(quantile(d$K27_DMSO, .9), 2), "\n")

# 2) positive controls: Fig 1C/1D boxed peaks (Eip74EF, Eip75B, br)
bp <- read.csv(inp(BOX_PEAKS))
bp <- bp[bp$gene %in% c("Eip74EF", "Eip75B", "br"), ]
y <- x[x$peak %in% bp$peak, ]
cat("\n== positive control: Fig 1C/D boxed peaks (n =", nrow(y), ")  H3K27ac lfc median", round(median(y$K27_lfc), 2),
    " up in both reps:", sum(y$K27_lfc_r1 > 0 & y$K27_lfc_r2 > 0), "  EcR lfc median", round(median(y$EcR_lfc), 2), "\n")

# 3) genome-wide: ATAC change vs H3K27ac change (distal peaks)
cat("\n== distal peaks: H3K27ac lfc by ATAC DAR group (1 h and 3 h)\n")
for (g in c("atac1", "atac3")) {
  t <- tapply(d$K27_lfc, d[[g]], function(v) c(n = length(v), median = round(median(v), 3)))
  print(do.call(rbind, t)[c("up", "unchanged", "down"), ])
  cat("  Wilcoxon up vs unchanged P =", signif(wilcox.test(d$K27_lfc[d[[g]] == "up"], d$K27_lfc[d[[g]] == "unchanged"])$p.value, 2),
      "; down vs unchanged P =", signif(wilcox.test(d$K27_lfc[d[[g]] == "down"], d$K27_lfc[d[[g]] == "unchanged"])$p.value, 2), "\n")
}
cat("Spearman ATAC 1 h log2FC vs H3K27ac lfc (distal):", round(cor(d$log2FoldChange_1h, d$K27_lfc, method = "s", use = "complete"), 3),
    "; ATAC 3 h:", round(cor(d$log2FoldChange_3h, d$K27_lfc, method = "s", use = "complete"), 3), "\n")
cat("replicate agreement of H3K27ac lfc (distal, Pearson):", round(cor(d$K27_lfc_r1, d$K27_lfc_r2), 3), "\n")
