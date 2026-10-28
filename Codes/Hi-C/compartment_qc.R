#!/usr/bin/env Rscript
# One-sided Wilcoxon QC of peak density and GC content in Hi-C compartments.
# One-sided Wilcoxon QC of Hi-C compartment peak density and GC content.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: Rscript compartment_qc.R <input_positive_counts.txt> <input_negative_counts.txt> <input_positive_gc.txt> <input_negative_gc.txt> <output_qc.tsv>")
}
read_column <- function(path, column) {
  if (!file.exists(path)) stop("Input not found: ", path)
  dat <- read.table(path, header = FALSE, sep = "\t", comment.char = "#")
  if (ncol(dat) < column) stop("Missing column ", column, " in ", path)
  x <- as.numeric(dat[[column]])
  x[is.finite(x)]
}
positive_peaks <- read_column(args[[1]], 4)
negative_peaks <- read_column(args[[2]], 4)
positive_gc <- read_column(args[[3]], 5)
negative_gc <- read_column(args[[4]], 5)
if (any(lengths(list(positive_peaks, negative_peaks, positive_gc, negative_gc)) == 0L)) {
  stop("One of the input groups is empty")
}
peak_test <- wilcox.test(positive_peaks, negative_peaks, alternative = "greater")
gc_test <- wilcox.test(positive_gc, negative_gc, alternative = "greater")
result <- data.frame(
  comparison = c("peak_density_positive_gt_negative", "gc_positive_gt_negative"),
  n_positive = c(length(positive_peaks), length(positive_gc)),
  n_negative = c(length(negative_peaks), length(negative_gc)),
  median_positive = c(median(positive_peaks), median(positive_gc)),
  median_negative = c(median(negative_peaks), median(negative_gc)),
  statistic_W = c(unname(peak_test$statistic), unname(gc_test$statistic)),
  p_value = c(peak_test$p.value, gc_test$p.value)
)
dir.create(dirname(args[[5]]), recursive = TRUE, showWarnings = FALSE)
write.table(result, args[[5]], sep = "\t", quote = FALSE, row.names = FALSE)
