#!/usr/bin/env Rscript
# Salmon/tximport and DESeq2 differential expression with apeglm shrinkage.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) stop("Usage: Rscript run_rna_deseq2.R input_samples.tsv input_tx2gene.tsv output_dir")
input_samples <- args[[1]]
input_tx2gene <- args[[2]]
output_dir <- args[[3]]
for (path in c(input_samples, input_tx2gene)) if (!file.exists(path)) stop("Input missing: ", path)
suppressPackageStartupMessages({library(tximport); library(DESeq2); library(apeglm); library(ggplot2)})
samples <- read.delim(input_samples, check.names=FALSE, stringsAsFactors=FALSE)
if (!all(c("sample_id", "path", "condition") %in% names(samples))) stop("Sample table needs sample_id, path, condition")
if (anyDuplicated(samples$sample_id) || !all(c("ctrl", "20E") %in% samples$condition)) stop("Invalid sample IDs or conditions")
quant_paths <- path.expand(samples$path)
relative <- !grepl("^/", quant_paths)
quant_paths[relative] <- file.path(dirname(normalizePath(input_samples)), quant_paths[relative])
if (!all(file.exists(quant_paths))) stop("Missing quant.sf: ", paste(quant_paths[!file.exists(quant_paths)], collapse=", "))
names(quant_paths) <- samples$sample_id
tx2gene <- read.delim(input_tx2gene, header=FALSE, col.names=c("tx_id", "gene_id"),
                      stringsAsFactors=FALSE)
txi <- tximport(files=quant_paths, type="salmon", tx2gene=tx2gene,
                ignoreTxVersion=TRUE)
col_data <- data.frame(condition=factor(samples$condition, levels=c("ctrl", "20E")),
                       row.names=samples$sample_id)
dds <- DESeqDataSetFromTximport(txi=txi, colData=col_data, design=~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
if (!"condition_20E_vs_ctrl" %in% resultsNames(dds)) stop("Expected 20E vs ctrl coefficient missing")
res <- lfcShrink(dds, coef="condition_20E_vs_ctrl", type="apeglm")
tab <- data.frame(gene=rownames(res), as.data.frame(res), row.names=NULL)
tab <- tab[order(tab$pvalue, na.last=TRUE), , drop=FALSE]
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
write.csv(tab, file.path(output_dir, "output_deseq2.csv"), row.names=FALSE)
png(file.path(output_dir, "output_dispersion.png"))
plotDispEsts(dds)
dev.off()
vsd <- vst(dds, blind=FALSE)
pca <- plotPCA(vsd, intgroup="condition") + ggtitle("20E vs control") + theme_bw()
ggsave(file.path(output_dir, "output_pca.png"), plot=pca, width=6, height=5)
cat("Genes:", nrow(tab), "Significant padj < 0.05:", sum(tab$padj < 0.05, na.rm=TRUE), "\n")
