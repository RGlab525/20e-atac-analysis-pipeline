# Fuzzy c-means of TSS DARs (0 h vs 3 h) into temporal groups T1-T4.
#
# Manuscript: Fig 3C, Supplementary Table S5 (T1 n = 68, T2 n = 63, T3 n = 825, T4 n = 443)
#   Up DARs and down DARs are clustered separately into 2 clusters each (e1071::cmeans, m = 2),
#   using log2FC at 0 h (set to 0), 1 h and 3 h.
#   Up cluster 1 = T1, up cluster 2 = T2, down cluster 1 = T3, down cluster 2 = T4 (with set.seed(123)).
# Input : input_wald_tss_3h.csv, input_wald_tss_1h.csv
# Output: $ATAC_OUT/02_cmeans/output_TSS_cmeans_groups.csv
#         (peak, X0h, X1h, X3h, direction, cluster, membership_1, membership_2, membership, group, gene)
# Source: 2025-03-10_ATAC_seq.Rmd lines 585-638 (DAR selection) and 1131-1588 (c-means, run of 2026-03-10).
# Changes from the original code:
#   - The original saved the membership matrix with `$membership <- cbind(...)`, which nested a copy of the
#     data frame into the CSV; here membership columns are written directly. Clusters are unchanged.
#   - Group labels T1-T4 are added (the original CSVs only had cluster 1/2).
source("config/config.R")
suppressMessages(library(e1071))
out <- out_dir("02_cmeans")

tss.3h <- read.delim2(inp("input_wald_tss_3h.csv"), header = TRUE, sep = ",", row.names = 1)
tss.1h <- read.delim2(inp("input_wald_tss_1h.csv"), header = TRUE, sep = ",", row.names = 1)

fc.3h <- as.numeric(tss.3h$log2FoldChange); padj.3h <- as.numeric(tss.3h$padj)
up.peaks   <- rownames(tss.3h)[which(fc.3h > log2(1.5)   & padj.3h < 0.05)]
down.peaks <- rownames(tss.3h)[which(fc.3h < log2(1/1.5) & padj.3h < 0.05)]

# rows in the order of the 1 h table, as in the original
data_matrix <- data.frame(
    X0h = 0,
    X1h = as.numeric(tss.1h[rownames(tss.1h), "log2FoldChange"]),
    X3h = as.numeric(tss.3h[rownames(tss.1h), "log2FoldChange"]),
    row.names = rownames(tss.1h))

run_cmeans <- function(m, labels, direction) {
    set.seed(123)
    cm <- cmeans(as.matrix(m), centers = 2, iter.max = 100, m = 2)
    data.frame(peak = rownames(m), m, direction = direction, cluster = cm$cluster,
               membership_1 = cm$membership[, 1], membership_2 = cm$membership[, 2],
               membership = cm$membership[cbind(seq_len(nrow(m)), cm$cluster)],
               group = labels[cm$cluster], row.names = NULL)
}

res <- rbind(
    run_cmeans(data_matrix[rownames(data_matrix) %in% up.peaks, ],   c("T1", "T2"), "up"),
    run_cmeans(data_matrix[rownames(data_matrix) %in% down.peaks, ], c("T3", "T4"), "down"))
res$gene <- tss.3h[res$peak, "gene"]

write.csv(res, file.path(out, "output_TSS_cmeans_groups.csv"), row.names = FALSE)
print(table(res$group))
