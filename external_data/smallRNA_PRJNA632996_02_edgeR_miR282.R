# =============================================================================
# R3 Major 8: response of mature miR-282-5p / -3p to 20E in the public S2 small RNA-seq (PRJNA632996)
#   each time point separately (6 h: DMSO 3 vs 20E 3; 2 d: DMSO 3 vs 20E 3), edgeR quasi-likelihood F test
#   filter: CPM >= 1 in >= 3 of the 6 samples (miR-282 always kept); normalization: miRNA-assigned reads + TMM
#   Used for the response letter and Supplementary Figure S7C.
# Input:  $ATAC_OUT/smallRNA_PRJNA632996/output_mature_miRNA_counts.tsv (step 01)
# Output: $ATAC_OUT/smallRNA_PRJNA632996/edgeR_<time>_all_miRNA.csv, output_miR282_summary.csv, output_miR282_CPM_per_sample.csv
# Source: 261016_revision/14_R3_Major8_miR282/scripts/02_edgeR_miR282.R (paths only changed)
# =============================================================================
source("config/config.R")
suppressMessages(library(edgeR))
O <- out_dir("smallRNA_PRJNA632996")
x <- read.delim(file.path(O, "output_mature_miRNA_counts.tsv"), row.names = 1, check.names = FALSE)
out <- list(); cpmL <- list()
for (tp in c("6h", "2d")) {
  m <- as.matrix(x[, grep(paste0("_", tp, "_"), colnames(x))])
  g <- factor(sub("_.*", "", colnames(m)), levels = c("DMSO", "20E"))
  y <- DGEList(m, group = g)
  keep <- rowSums(cpm(y) >= 1) >= 3 | grepl("miR-282", rownames(y))
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  d <- model.matrix(~ g); y <- estimateDisp(y, d); f <- glmQLFit(y, d)
  t <- topTags(glmQLFTest(f, coef = 2), n = Inf, sort.by = "PValue")$table
  t$rank_by_P <- seq_len(nrow(t)); t$n_tested <- nrow(t)
  write.csv(t, file.path(O, sprintf("output_edgeR_%s_all_miRNA.csv", tp)))
  for (n in c("dme-miR-282-5p", "dme-miR-282-3p")) out[[paste(tp, n)]] <- data.frame(time = tp, miRNA = n, t[n, c("logFC", "logCPM", "PValue", "FDR", "rank_by_P", "n_tested")])
  c1 <- cpm(y); cpmL[[tp]] <- data.frame(time = tp, sample = colnames(c1), group = g, t(c1[c("dme-miR-282-5p", "dme-miR-282-3p"), ]), check.names = FALSE)
}
s <- do.call(rbind, out); rownames(s) <- NULL
write.csv(s, file.path(O, "output_miR282_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, cpmL), file.path(O, "output_miR282_CPM_per_sample.csv"), row.names = FALSE)
print(s); print(do.call(rbind, cpmL))
