# =============================================================================
# Fig 5F (revised) and Supplementary Fig S5E: RNA-seq (GSE47691, 24 h) log2FC of genes in the ATAC c-means groups
#   ATAC part of 5F: peak log2FC at 1 h and 3 h of T1 / T2 (Table S5), medians only (no test)
#   RNA: Gene column of Tables S5 / S8 -> FlyBase ID (org.Dm.eg.db SYMBOL, then ALIAS for unmatched names),
#        genes present in the RNA-seq DESeq2 table, each gene counted once per group
#   test: two-sided Wilcoxon rank-sum, T1 vs T2, T3 vs T4, E1 vs E2, E3 vs E4; BH over the 4 comparisons;
#         for the enhancer comparisons genes shared by the two groups are removed
# Input:  supplementary_tables/input_tss_groups.xlsx, input_enhancer_groups.xlsx,
#         input_rna24h_deseq2.csv
# Output: $ATAC_OUT/Fig5F/output_Fig5F_RNA_24h_wilcoxon.csv, output_Fig5F_ATAC_T1T2_medians.csv
# Source: 261016_revision/06_Fig5_corrections/scripts/fig5F_T1T2_ATAC_RNA.R (plotting left out)
# =============================================================================
source("config/config.R")
suppressMessages({ library(readxl); library(org.Dm.eg.db) })
out <- out_dir("Fig5F")

rna <- read.csv(inp("input_rna24h_deseq2.csv"), row.names = 1); rna <- rna[!is.na(rna$log2FoldChange), ]
rd <- function(f) { d <- as.data.frame(suppressMessages(read_excel(inp("supplementary_tables", f), skip = 1))); d[!is.na(d$Cluster), ] }
s5 <- rd("input_tss_groups.xlsx"); cl <- rbind(s5, rd("input_enhancer_groups.xlsx"))
sp <- function(x) { g <- trimws(unlist(strsplit(as.character(x), "[;,]"))); unique(g[!is.na(g) & !g %in% c("", "NA", "-")]) }
al <- unique(unlist(lapply(cl$Gene, sp)))
fb <- suppressMessages(mapIds(org.Dm.eg.db, al, "FLYBASE", "SYMBOL", multiVals = "first"))
m <- is.na(fb); fb[m] <- suppressMessages(mapIds(org.Dm.eg.db, al[m], "FLYBASE", "ALIAS", multiVals = "first"))
G <- lapply(setNames(nm = c(paste0("T", 1:4), paste0("E", 1:4))), function(g) {
    f <- unique(na.omit(fb[sp(cl$Gene[cl$Cluster == g])])); f[f %in% rownames(rna)] })

st <- do.call(rbind, lapply(list(c("T1", "T2"), c("T3", "T4"), c("E1", "E2"), c("E3", "E4")), function(p) {
    a <- G[[p[1]]]; b <- G[[p[2]]]
    if (substr(p[1], 1, 1) == "E") { sh <- intersect(a, b); a <- setdiff(a, sh); b <- setdiff(b, sh) }
    x <- rna[a, "log2FoldChange"]; y <- rna[b, "log2FoldChange"]
    data.frame(comparison = paste(p, collapse = " vs "), n1 = length(x), n2 = length(y), median1 = median(x), median2 = median(y),
               p = wilcox.test(x, y)$p.value) }))
st$q_BH <- p.adjust(st$p, "BH")
write.csv(st, file.path(out, "output_Fig5F_RNA_24h_wilcoxon.csv"), row.names = FALSE); print(st, digits = 3)

t12 <- s5[s5$Cluster %in% c("T1", "T2"), ]
med <- aggregate(cbind(`Log2FC (1h)`, `Log2FC (3h)`) ~ Cluster, t12, median)
med$n_peaks <- as.integer(table(t12$Cluster)[med$Cluster])
write.csv(med, file.path(out, "output_Fig5F_ATAC_T1T2_medians.csv"), row.names = FALSE); print(med)
