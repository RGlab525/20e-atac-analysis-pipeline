# =============================================================================
# Fig 5C, 5D (revised): temporal ATAC c-means groups (T1-T4, E1-E4) vs RNA-seq DEGs (GSE47691, 24 h)
#   one-sided Fisher's exact test (greater) of each group's genes among RNA up / down genes;
#   BH correction within TSS groups (8 tests) and within enhancer groups (8 tests)
#   unit = unique gene (a gene with several peaks in a group is counted once)
#   RNA DEG: padj < 0.05 & |log2FC| > log2(1.5); background = all RNA-seq genes with a symbol
#   Group genes = Gene column of Supplementary Tables S5 (T1-T4 = 68/63/825/443 peaks) and S8
# Input:  supplementary_tables/input_tss_groups.xlsx, input_enhancer_groups.xlsx,
#         input_rna24h_deseq2.csv (RNA-seq DESeq2)
# Output: $ATAC_OUT/Fig5CD/output_Fig5CD_fisher.csv
# Source: 261016_revision/06_Fig5_corrections/scripts/fig5_CFG_revision.R, sections 1-3 (plotting left out)
# =============================================================================
source("config/config.R")
suppressMessages({ library(readxl); library(org.Dm.eg.db) })
outdir <- out_dir("Fig5CD")

# ----- 1. RNA-seq -----
rna <- read.csv(inp("input_rna24h_deseq2.csv"),
                stringsAsFactors = FALSE, row.names = 1)
rna$symbol <- suppressMessages(mapIds(org.Dm.eg.db, keys = rownames(rna),
                     column = "SYMBOL", keytype = "FLYBASE", multiVals = "first"))
rna <- rna[!is.na(rna$symbol) & !duplicated(rna$symbol), ]

fc.cut   <- log2(1.5)
padj.cut <- 0.05
rna.bg   <- rna$symbol
rna.up   <- rna$symbol[which(rna$padj < padj.cut & rna$log2FoldChange >  fc.cut)]
rna.down <- rna$symbol[which(rna$padj < padj.cut & rna$log2FoldChange < -fc.cut)]
cat(sprintf("RNA-seq detected: %d  up: %d  down: %d\n",
            length(rna.bg), length(rna.up), length(rna.down)))

# ----- 2. ATAC cluster -> gene (Gene column of Supplementary Tables S5 / S8) -----
read.supp <- function(f) as.data.frame(suppressMessages(read_excel(inp("supplementary_tables", f), skip = 1)))
s5 <- read.supp("input_tss_groups.xlsx")
s8 <- read.supp("input_enhancer_groups.xlsx")

split.genes <- function(x) {
    g <- trimws(unlist(strsplit(as.character(x), "[;,]")))
    unique(g[!is.na(g) & g != "" & g != "NA"])
}
groups <- c("T1", "T2", "T3", "T4", "E1", "E2", "E3", "E4")
group.peaks <- lapply(setNames(groups, groups), function(g) {
    d <- if (substr(g, 1, 1) == "T") s5 else s8
    d$`Peak ID`[d$Cluster == g]
})
group.genes <- lapply(setNames(groups, groups), function(g) {
    d <- if (substr(g, 1, 1) == "T") s5 else s8
    intersect(split.genes(d$Gene[d$Cluster == g]), rna.bg)
})

# ----- 3. Fisher (Fig 5C / 5D) -----
run.fisher <- function(grp.g, dir.g, bg) {
    a <- length(intersect(grp.g, dir.g))
    b <- length(grp.g) - a
    c <- length(intersect(setdiff(bg, grp.g), dir.g))
    d <- length(bg) - a - b - c
    ft <- fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE), alternative = "greater")
    c(n_overlap = a, n_grp_only = b, n_rna_only = c, n_neither = d,
      odds_ratio = unname(ft$estimate), pvalue = ft$p.value)
}
fisher <- do.call(rbind, lapply(groups, function(g) do.call(rbind, lapply(c("up", "down"), function(dn) {
    r <- run.fisher(group.genes[[g]], if (dn == "up") rna.up else rna.down, rna.bg)
    data.frame(group = g, rna_dir = dn, n_peaks = length(group.peaks[[g]]),
               n_genes = length(group.genes[[g]]), t(r))
}))))
fisher$set <- ifelse(substr(fisher$group, 1, 1) == "T", "TSS", "ENH")
fisher$padj_BH <- ave(fisher$pvalue, fisher$set, FUN = function(p) p.adjust(p, "BH"))
write.csv(fisher, file.path(outdir, "output_Fig5CD_fisher.csv"), row.names = FALSE)
print(transform(fisher[, c("group", "rna_dir", "n_genes", "n_overlap", "odds_ratio", "padj_BH")], odds_ratio = round(odds_ratio, 2), padj_BH = signif(padj_BH, 3)), row.names = FALSE)
