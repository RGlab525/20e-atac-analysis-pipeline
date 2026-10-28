# Fuzzy c-means of enhancer DARs (non-TSS peaks, 0 h vs 3 h) into temporal groups E1-E4.
#
# Manuscript: Fig 4C, Supplementary Table S8 (E1 n = 486, E2 n = 792, E3 n = 456, E4 n = 178)
#   Same procedure as the TSS groups. With set.seed(123):
#   up cluster 2 = E1, up cluster 1 = E2, down cluster 1 = E3, down cluster 2 = E4.
# Input : input_peak_annotations.xlsx, sheet "260322_1h_3h"
#         (log2FC / padj at 1 h and 3 h in one row per peak; TSS = "X" marks non-TSS peaks)
# Output: $ATAC_OUT/02_cmeans/output_enhancer_cmeans_groups.csv
#         (peak, X0h, X1h, X3h, direction, cluster, membership_1, membership_2, membership, group, linked_genes)
# Source: 2025-03-10_ATAC_seq.Rmd lines 1592-2024 (run of 2026-07-01, files 260701_cmean_cluster_ENH_*.csv).
#         The E1-E4 labels were added to those CSVs by hand (E_group column); here they are added in code.
source("config/config.R")
suppressMessages({ library(e1071); library(readxl) })
out <- out_dir("02_cmeans")

enh.tab <- as.data.frame(read_excel(inp("input_peak_annotations.xlsx"), sheet = "260322_1h_3h"))
rownames(enh.tab) <- enh.tab$peak_id
for (k in c("log2FoldChange_1h", "log2FoldChange_3h", "padj_3h")) enh.tab[[k]] <- as.numeric(enh.tab[[k]])
enh.tab <- enh.tab[enh.tab$TSS == "X", ]

up   <- enh.tab[enh.tab$log2FoldChange_3h >  log2(1.5)   & enh.tab$padj_3h < 0.05, ]
down <- enh.tab[enh.tab$log2FoldChange_3h <  log2(1/1.5) & enh.tab$padj_3h < 0.05, ]

data_matrix <- data.frame(X0h = 0, X1h = enh.tab$log2FoldChange_1h, X3h = enh.tab$log2FoldChange_3h,
                          row.names = rownames(enh.tab))

run_cmeans <- function(m, labels, direction) {
    set.seed(123)
    cm <- cmeans(as.matrix(m), centers = 2, iter.max = 100, m = 2)
    data.frame(peak = rownames(m), m, direction = direction, cluster = cm$cluster,
               membership_1 = cm$membership[, 1], membership_2 = cm$membership[, 2],
               membership = cm$membership[cbind(seq_len(nrow(m)), cm$cluster)],
               group = labels[cm$cluster], row.names = NULL)
}

res <- rbind(
    run_cmeans(na.omit(data_matrix[rownames(up), ]),   c("E2", "E1"), "up"),
    run_cmeans(na.omit(data_matrix[rownames(down), ]), c("E3", "E4"), "down"))
res$linked_genes <- enh.tab[res$peak, "enhancer"]

write.csv(res, file.path(out, "output_enhancer_cmeans_groups.csv"), row.names = FALSE)
print(table(res$group))
