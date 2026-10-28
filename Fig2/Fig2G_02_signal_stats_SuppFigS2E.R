# =============================================================================
# Fig 2G and Supplementary Fig. S2E: statistics of the mean ATAC signal per locus over time
#   Signal = RPGC bigWig mean per window (multiBigwigSummary, Fig2AB_02) or per region (computeMatrix, SuppFig_S2_03)
#   Statistics (no multiple-testing correction; one test per time point):
#     - each time point vs 0 h: paired Wilcoxon signed-rank test across loci
#     - TSS peaks vs distal peaks: per-locus log2((t + 0.1) / (0 h + 0.1)) compared by Wilcoxon rank-sum test
#   Changes are shown as log2 ratios because TSS peaks start from a higher 0-h signal (median 5.1 vs 3.1).
#   Unit = locus; the two replicates are technical, so these are not tests across biological replicates.
#   Windows without bigWig signal are set to 0.
#
# Input:
#   $ATAC_OUT/Fig2AB/boxplot/{TSS_peaks,ENH_distal9229}_mean_signal.tab
#   $ATAC_OUT/SuppFig_S2/boxplot/{Intergenic,Intergenic_noPeak,GeneBody,GeneBody_noPeak}_mean_signal.tab
# Output ($ATAC_OUT/Fig2G/):
#   output_Fig2_boxplot_stats_vs_0h.csv, output_Fig2_TSSpeaks_vs_distalPeaks_comparison.csv  (numbers in Results and the response letter)
#
# Source: 261016_revision/03_Fig2AB_heatmap_boxplot/scripts/03_boxplot_stats.R
#   (box plots and the original-set comparison panels are left out)
# =============================================================================
source("config/config.R")
O   <- out_dir("Fig2G")
tp   <- TIMES
sets <- c(TSS_peaks = "TSS peaks\n(new Fig 2A)", ENH_distal9229 = "Distal peaks\n(new Fig 2B)",
          Intergenic = "Intergenic\n(same regions and values as Supp. Fig. S2A)",
          Intergenic_noPeak = "Intergenic\nwithout ATAC peaks",
          GeneBody = "Gene body\n(TSS +1 kb to TES)", GeneBody_noPeak = "Gene body\nwithout ATAC peaks")
src <- c(TSS_peaks = "Fig2AB/boxplot", ENH_distal9229 = "Fig2AB/boxplot", Intergenic = "SuppFig_S2/boxplot",
         Intergenic_noPeak = "SuppFig_S2/boxplot", GeneBody = "SuppFig_S2/boxplot", GeneBody_noPeak = "SuppFig_S2/boxplot")

rd <- function(s) {
    x <- read.delim(file.path(OUT_ROOT, src[[s]], paste0("output_", s, "_mean_signal.tab")), check.names = FALSE)
    m <- as.matrix(x[, 4:9]); colnames(m) <- tp; m[is.na(m)] <- 0; m
}
ps <- 0.1; lfc <- function(m) log2((m + ps) / (m[, "0h"] + ps))     # per-locus log2(t / 0 h)
star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "n.s.")))
vals <- lapply(setNames(names(sets), names(sets)), rd)

st <- do.call(rbind, lapply(names(vals), function(s) { m <- vals[[s]]
    do.call(rbind, lapply(tp, function(t) { d <- lfc(m)[, t]
        p <- if (t == "0h") NA else wilcox.test(m[, t], m[, "0h"], paired = TRUE)$p.value
        data.frame(set = s, n = nrow(m), time = t, median = median(m[, t]), IQR = IQR(m[, t]),
                   median_log2FC_vs_0h = median(d), pct_change_of_median = 100 * (median(m[, t]) / median(m[, "0h"]) - 1),
                   pct_loci_decreased = if (t == "0h") NA else 100 * mean(d < 0), p_paired_wilcoxon_vs_0h = p,
                   star = if (t == "0h") "" else star(p)) })) }))
write.csv(st, file.path(O, "output_Fig2_boxplot_stats_vs_0h.csv"), row.names = FALSE)

cmp <- do.call(rbind, lapply(tp[-1], function(t) {
    a <- lfc(vals$TSS_peaks)[, t]; b <- lfc(vals$ENH_distal9229)[, t]
    data.frame(time = t, TSS_median_log2FC = median(a), ENH_median_log2FC = median(b),
               difference_of_medians = median(a) - median(b), p_rank_sum_log2FC = wilcox.test(a, b)$p.value) }))
cmp$star <- star(cmp$p_rank_sum_log2FC)
cmp$p_text <- ifelse(cmp$p_rank_sum_log2FC == 0, "< 1e-300 (below the smallest number R reports)", format(signif(cmp$p_rank_sum_log2FC, 3)))
write.csv(cmp, file.path(O, "output_Fig2_TSSpeaks_vs_distalPeaks_comparison.csv"), row.names = FALSE)
print(transform(st, median = round(median, 3), median_log2FC_vs_0h = round(median_log2FC_vs_0h, 3),
                pct_change_of_median = round(pct_change_of_median, 1), pct_loci_decreased = round(pct_loci_decreased, 1),
                p_paired_wilcoxon_vs_0h = signif(p_paired_wilcoxon_vs_0h, 3), IQR = round(IQR, 3)), row.names = FALSE)
print(transform(cmp, p_rank_sum_log2FC = signif(p_rank_sum_log2FC, 3)), row.names = FALSE)
