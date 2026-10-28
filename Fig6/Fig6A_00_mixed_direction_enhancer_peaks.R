# =============================================================================
# Fig 6A / Supplementary Fig S6, step 0: genes linked to both significantly increased and significantly decreased
#   enhancer (distal, TSS == "X") peaks at 3 h; all distal peaks linked to these genes are listed
#   link = "enhancer" column of the peak table (10-kb assignment); direction at 3 h: padj_3h < 0.05 and
#   log2FC_3h > log2(1.5) (up) or < -log2(1.5) (down), otherwise ns
#
# Input:  PEAK_TABLE, sheet "260322_1h_3h" (peak_id, chr, start, end, log2FoldChange_1h, log2FoldChange_3h, padj_3h,
#         TSS, enhancer). Default = input_fig6_peak_annotations.xlsx, the table the published Fig 6A and S6
#         were made from.
# Output: $ATAC_OUT/Fig6A/output_mixed_direction_peaks.csv
#         (peak_id, chr, start, end, log2FC_1h, log2FC_3h, padj_3h, gene, direction_3h; log2FC rounded to 4 decimals)
#
# Source: written for the code release; the analysed file output_mixed_direction_peaks.csv was made outside the
#   Rmd. With the default table this script reproduces it exactly (311 peak-gene rows, 72 genes).
# Note: in the 260322 table padj_3h is the BH adjustment over all 14,003 peaks (1,095 up / 589 down distal DARs).
#   The final table used for Fig 4 (input_peak_annotations.xlsx; 1,278 up / 634 down) has
#   different padj values for distal peaks (same log2FC); with PEAK_TABLE set to that file the list becomes
#   382 rows / 91 genes. Fig 6A / S6 have not been regenerated with the final table.
# =============================================================================
source("config/config.R")
suppressMessages(library(readxl))
PEAK_TABLE <- Sys.getenv("PEAK_TABLE", "input_fig6_peak_annotations.xlsx")
O <- out_dir("Fig6A")

a <- suppressWarnings(as.data.frame(read_excel(inp(PEAK_TABLE), sheet = "260322_1h_3h")))
d <- a[a$TSS == "X", ]
L <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    g <- trimws(strsplit(as.character(d$enhancer[i]), "[;,]")[[1]]); g <- g[!g %in% c("-", "", NA)]
    if (!length(g)) return(NULL)
    data.frame(peak_id = d$peak_id[i], chr = d$chr[i], start = d$start[i], end = d$end[i],
               log2FC_1h = round(as.numeric(d$log2FoldChange_1h[i]), 4), log2FC_3h = round(as.numeric(d$log2FoldChange_3h[i]), 4),
               padj_3h = as.numeric(d$padj_3h[i]), gene = g)
}))
L$direction_3h <- ifelse(!is.na(L$padj_3h) & L$padj_3h < PADJ_CUT & L$log2FC_3h >  FC_CUT, "up",
                  ifelse(!is.na(L$padj_3h) & L$padj_3h < PADJ_CUT & L$log2FC_3h < -FC_CUT, "down", "ns"))
mixed <- names(which(tapply(L$direction_3h, L$gene, function(v) any(v == "up") & any(v == "down"))))
L <- L[L$gene %in% mixed, ]; L <- L[order(L$gene), ]
write.csv(L, file.path(O, "output_mixed_direction_peaks.csv"), row.names = FALSE)
cat("peak-gene rows:", nrow(L), " genes:", length(mixed), "\n"); print(table(L$direction_3h))
