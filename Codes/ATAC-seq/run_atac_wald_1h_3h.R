#!/usr/bin/env Rscript
# DESeq2 Wald contrasts: 1 h versus 0 h and 3 h versus 0 h.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4L) stop("Usage: Rscript run_atac_wald_1h_3h.R input_counts.csv input_metadata.tsv input_tss_peaks.csv output_dir")
input_counts <- args[[1]]
input_metadata <- args[[2]]
input_tss <- args[[3]]
output_dir <- args[[4]]
for (path in c(input_counts, input_metadata, input_tss)) if (!file.exists(path)) stop("Input missing: ", path)
suppressPackageStartupMessages(library(DESeq2))
counts <- read.csv(input_counts, row.names=1, check.names=FALSE)
metadata <- read.table(input_metadata, header=TRUE, row.names=1, check.names=FALSE,
                       stringsAsFactors=FALSE)
tss <- read.csv(input_tss, stringsAsFactors=FALSE, check.names=FALSE)
if (!"peak_id" %in% names(tss)) stop("TSS table needs peak_id column")
if (!all(c("time", "batch") %in% names(metadata))) stop("Metadata needs time and batch columns")
metadata <- metadata[metadata$time %in% c("0h", "1h", "3h"), , drop=FALSE]
if (!all(rownames(metadata) %in% colnames(counts))) stop("Metadata samples missing from count matrix")
counts <- counts[, rownames(metadata), drop=FALSE]
if (anyNA(counts) || any(counts < 0) || any(counts != floor(counts))) stop("Counts must be nonnegative integers")
if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) stop("Duplicate peak or sample IDs")
metadata$time <- factor(metadata$time, levels=c("0h", "1h", "3h"))
metadata$batch <- factor(metadata$batch)
if (any(table(metadata$time) == 0L) || nlevels(metadata$batch) < 2L) stop("Missing time group or batch")
dds <- DESeqDataSetFromMatrix(countData=as.matrix(counts), colData=metadata,
                              design=~batch+time)
dds <- DESeq(dds, test="Wald")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
tss_ids <- unique(tss$peak_id)
for (hour in c("1h", "3h")) {
  res <- results(dds, contrast=c("time", hour, "0h"))
  tab <- data.frame(peak_id=rownames(res), as.data.frame(res), row.names=NULL,
                    check.names=FALSE)
  write.csv(tab, file.path(output_dir, paste0("output_wald_", hour, "_all.csv")), row.names=FALSE)
  tss_tab <- tab[tab$peak_id %in% tss_ids, , drop=FALSE]
  if ("matched_gene" %in% names(tss)) {
    tss_tab$matched_gene <- tss$matched_gene[match(tss_tab$peak_id, tss$peak_id)]
  }
  enh_tab <- tab[!tab$peak_id %in% tss_ids, , drop=FALSE]
  write.csv(tss_tab, file.path(output_dir, paste0("output_wald_", hour, "_tss.csv")), row.names=FALSE)
  write.csv(enh_tab, file.path(output_dir, paste0("output_wald_", hour, "_enhancer.csv")), row.names=FALSE)
  cat(hour, "all", nrow(tab), "TSS", nrow(tss_tab), "enhancer", nrow(enh_tab), "\n")
}
