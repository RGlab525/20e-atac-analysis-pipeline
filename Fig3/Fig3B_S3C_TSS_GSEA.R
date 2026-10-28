# Gene set enrichment analysis (GO Biological Process) of genes ranked by TSS accessibility change, 0 h vs 3 h.
#
# Manuscript: Fig 3B (terms with NES > 0), Supplementary Fig S3C (NES < 0),
#             Supplementary Tables S3 (30 terms, NES > 0) and S4 (47 terms, NES < 0); adjusted P < 0.05
# Ranking: DESeq2 Wald statistic (log2FC / lfcSE) of each TSS peak; a peak assigned to several genes
#          ("gene1; gene2") is counted for each gene; for a gene with several peaks the largest |stat| is kept.
# Input : input_wald_tss_3h.csv
# Output: $ATAC_OUT/Fig3B/output_TSS_gseGO_BP_all.tsv              all tested terms (core genes as symbols)
#         $ATAC_OUT/Fig3B/TSS_gseGO_BP_padj0.05_{up,down}.tsv  = Tables S3 / S4 and the ClueGO input
# Network: the significant terms were loaded into Cytoscape (v3.9.1) / ClueGO (v2.5.10) to draw the
#          networks in Fig 3B and Supplementary Fig S3C (interactive; no script).
# Source: 2025-03-10_ATAC_seq.Rmd lines 3699-3836 (run of 2026-03-09, file 2026-03-09_BP_deseq.0hvs3h_newDMEdataset_symbol.txt).
#   Note: the analysis uses clusterProfiler::gseGO (which runs fgsea internally) with org.Dm.eg.db GO annotation,
#   not msigdbr gene sets. The msigdbr/fgsea chunk at Rmd lines 3435-3697 is an earlier version that was not used.
source("config/config.R")
suppressMessages({ library(clusterProfiler); library(org.Dm.eg.db) })
out <- out_dir("Fig3B")

tss <- read.delim2(inp("input_wald_tss_3h.csv"), header = TRUE, sep = ",", row.names = 1)

# one row per (peak, gene)
expanded <- do.call(rbind, lapply(seq_len(nrow(tss)), function(i) {
    genes <- trimws(strsplit(as.character(tss$gene[i]), "; ")[[1]])
    data.frame(gene = genes, stat = as.numeric(tss$stat[i]), stringsAsFactors = FALSE)
}))
expanded <- expanded[!is.na(expanded$gene) & expanded$gene != "", ]
expanded <- expanded[order(-abs(expanded$stat)), ]
expanded <- expanded[!duplicated(expanded$gene), ]
ranked.symbol <- sort(setNames(expanded$stat, expanded$gene), decreasing = TRUE)

gene.df <- bitr(names(ranked.symbol), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db)
gene.df <- gene.df[!is.na(gene.df$ENTREZID), ]
gene.df <- gene.df[!duplicated(gene.df$ENTREZID), ]
ranked <- ranked.symbol[gene.df$SYMBOL]
names(ranked) <- gene.df$ENTREZID
ranked <- sort(ranked[!is.na(ranked)], decreasing = TRUE)

set.seed(42)
gse <- gseGO(geneList = ranked, OrgDb = org.Dm.eg.db, keyType = "ENTREZID", ont = "BP",
             minGSSize = 3, maxGSSize = 1000, pvalueCutoff = 1, verbose = FALSE)

res <- as.data.frame(gse)
res <- res[order(res$p.adjust), ]
id2sym <- setNames(gene.df$SYMBOL, gene.df$ENTREZID)
res$core_enrichment <- vapply(res$core_enrichment, function(x) paste(id2sym[strsplit(x, "/")[[1]]], collapse = "/"), character(1))

write.table(res, file.path(out, "output_TSS_gseGO_BP_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sig <- res[res$p.adjust < 0.05, ]
write.table(sig[sig$NES > 0, ], file.path(out, "output_TSS_gseGO_BP_padj0.05_up.tsv"),   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig[sig$NES < 0, ], file.path(out, "output_TSS_gseGO_BP_padj0.05_down.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("significant terms: NES > 0", sum(sig$NES > 0), "; NES < 0", sum(sig$NES < 0), "\n")
