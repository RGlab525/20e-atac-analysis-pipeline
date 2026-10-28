# =============================================================================
# R1 Major 2 (1 h RNA-seq), step 3: 1 h RNA response (GSE304853) and comparison with our ATAC-seq and the 24 h RNA-seq
#   gene values from step 2 are median-of-ratios normalized between samples (size factors from genes with coverage >= 0.5
#   in all four samples); expressed = coverage >= 0.2 in any sample
#   1 h responsive gene (approximation, not a statistical test): both replicates change >= 1.5-fold in the same direction
#   comparison: TSS-peak genes, ATAC (DESeq2 1 h / 3 h) vs RNA 1 h (this step) vs RNA 24 h (GSE47691)
# Input:  $ATAC_OUT/RNA_1h_GSE304853/output_RNA_1h_gene_exon_coverage.csv (step 2), input_rna24h_deseq2.csv,
#         input_peak_annotations.xlsx (sheet 2)
# Output: $ATAC_OUT/RNA_1h_GSE304853/output_RNA_1h_gene_log2FC.csv, output_TSS_peak_genes_ATAC_vs_RNA_1h_24h.csv; statistics printed
#         (values shown in Supplementary Fig S8)
# Source: 261016_revision/13_R1_Major2_RNA_1h/scripts/03_RNA_1h_vs_ATAC.R (paths only changed)
# =============================================================================
source("config/config.R")
suppressMessages({ library(readxl); library(org.Dm.eg.db) }); options(warn = -1)
O <- out_dir("RNA_1h_GSE304853")
cv <- read.csv(file.path(O, "output_RNA_1h_gene_exon_coverage.csv"))
v <- cv[, c("DMSO_1", "DMSO_2", "E20_1", "E20_2")]
# between-sample normalization: median-of-ratios (size factors from genes with coverage >= 0.5 in all four samples)
ok <- apply(v, 1, function(z) all(z >= 0.5)); lg <- log(as.matrix(v[ok, ]))
sf <- exp(apply(lg - rowMeans(lg), 2, median)); v <- sweep(v, 2, sf, "/")
cat("size factors:", round(sf, 3), " (from", sum(ok), "genes)\n")
keep <- apply(v, 1, max) >= 0.2                                             # expressed genes (coverage >= 0.2 in any sample)
r1 <- data.frame(fbgn = cv$fbgn, symbol = cv$symbol, meanDMSO = (v$DMSO_1 + v$DMSO_2) / 2, mean20E = (v$E20_1 + v$E20_2) / 2,
                 lfc_r1 = log2((v$E20_1 + .1) / (v$DMSO_1 + .1)), lfc_r2 = log2((v$E20_2 + .1) / (v$DMSO_2 + .1)))[keep, ]
r1$log2FoldChange <- log2((r1$mean20E + .1) / (r1$meanDMSO + .1))
fc <- log2(1.5)
r1$DE1h <- ifelse(r1$lfc_r1 > fc & r1$lfc_r2 > fc, "up", ifelse(r1$lfc_r1 < -fc & r1$lfc_r2 < -fc, "down", "ns"))
write.csv(r1, file.path(O, "output_RNA_1h_gene_log2FC.csv"), row.names = FALSE)
cat("expressed genes:", nrow(r1), "  1 h DE up:", sum(r1$DE1h == "up"), " down:", sum(r1$DE1h == "down"), "\n")

# positive controls: early ecdysone genes
ctrl <- c("Eip74EF", "Eip75B", "br", "Hr3", "Hr4", "Hr39", "E23", "Kr-h1", "Blimp-1", "ImpL2", "Eip78C", "ftz-f1", "Eip93F", "lncRNA:CR43626")
z <- r1[match(ctrl, r1$symbol), c("symbol", "meanDMSO", "mean20E", "lfc_r1", "lfc_r2", "DE1h")]; z$symbol <- ctrl; z[, 2:5] <- round(z[, 2:5], 2); print(z, row.names = FALSE)

# 24 h RNA (manuscript)
r24 <- read.csv(inp("input_rna24h_deseq2.csv"), row.names = 1); r24$fbgn <- rownames(r24)
r24$DE24 <- ifelse(!is.na(r24$padj) & r24$padj < 0.05 & r24$log2FoldChange > fc, "up", ifelse(!is.na(r24$padj) & r24$padj < 0.05 & r24$log2FoldChange < -fc, "down", "ns"))

# ATAC: TSS peak -> gene (matched_gene symbol); symbol -> FBgn
a <- as.data.frame(read_excel(inp("input_peak_annotations.xlsx"), sheet = 2))
sy <- toTable(org.Dm.egSYMBOL); fbt <- toTable(org.Dm.egFLYBASE)
s2f <- function(s) fbt$flybase_id[match(sy$gene_id[match(s, sy$symbol)], fbt$gene_id)]
t <- a[a$TSS == "O" & a$matched_gene != "-", ]
t <- do.call(rbind, lapply(seq_len(nrow(t)), function(i) { g <- trimws(strsplit(t$matched_gene[i], "[;,]")[[1]])
  data.frame(gene = g, atac1 = t$log2FoldChange_1h[i], padj1 = t$padj_1h[i], atac3 = t$log2FoldChange_3h[i], padj3 = t$padj_3h[i]) }))
t <- t[!duplicated(t$gene), ]; t$fbgn <- s2f(t$gene)
t$ATAC3 <- ifelse(t$padj3 < 0.05 & t$atac3 > fc, "up", ifelse(t$padj3 < 0.05 & t$atac3 < -fc, "down", "ns")); t$ATAC3[is.na(t$ATAC3)] <- "ns"
m <- merge(t, r1[, c("fbgn", "log2FoldChange", "DE1h")], by = "fbgn"); names(m)[names(m) == "log2FoldChange"] <- "rna1"
m <- merge(m, r24[, c("fbgn", "log2FoldChange", "DE24")], by = "fbgn", all.x = TRUE); names(m)[names(m) == "log2FoldChange"] <- "rna24"
write.csv(m, file.path(O, "output_TSS_peak_genes_ATAC_vs_RNA_1h_24h.csv"), row.names = FALSE)
cat("\nTSS-peak genes with RNA data:", nrow(m), "(24 h available:", sum(!is.na(m$rna24)), ")\n")
cat(sprintf("Spearman ATAC (TSS peak) vs RNA:  1 h ATAC ~ 1 h RNA %.3f | 3 h ATAC ~ 1 h RNA %.3f | 3 h ATAC ~ 24 h RNA %.3f\n",
            cor(m$atac1, m$rna1, method = "s", use = "complete"), cor(m$atac3, m$rna1, method = "s", use = "complete"), cor(m$atac3, m$rna24, method = "s", use = "complete")))
cat("\nmedian RNA log2FC by 3 h TSS-peak DAR group (1 h RNA | 24 h RNA)\n")
for (g in c("up", "ns", "down")) { z <- m[m$ATAC3 == g, ]
  cat(sprintf("  ATAC %-4s n = %4d   1 h %+.3f   24 h %+.3f\n", g, nrow(z), median(z$rna1, na.rm = TRUE), median(z$rna24, na.rm = TRUE))) }
w <- function(col) c(up = signif(wilcox.test(m[[col]][m$ATAC3 == "up"], m[[col]][m$ATAC3 == "ns"])$p.value, 2),
                     down = signif(wilcox.test(m[[col]][m$ATAC3 == "down"], m[[col]][m$ATAC3 == "ns"])$p.value, 2))
cat("  Wilcoxon vs ns, 1 h RNA:", w("rna1"), "  24 h RNA:", w("rna24"), "\n")
cat("\nATAC 3 h up TSS genes that are RNA-up: 1 h", sum(m$ATAC3 == "up" & m$DE1h == "up"), "/", sum(m$ATAC3 == "up"),
    "; 24 h", sum(m$ATAC3 == "up" & m$DE24 == "up", na.rm = TRUE), "/", sum(m$ATAC3 == "up" & !is.na(m$DE24)), "\n")
cat("ATAC 3 h down TSS genes that are RNA-down: 1 h", sum(m$ATAC3 == "down" & m$DE1h == "down"), "/", sum(m$ATAC3 == "down"),
    "; 24 h", sum(m$ATAC3 == "down" & m$DE24 == "down", na.rm = TRUE), "/", sum(m$ATAC3 == "down" & !is.na(m$DE24)), "\n")
# agreement between 1 h and 24 h RNA
b <- merge(r1[, c("fbgn", "log2FoldChange")], r24[, c("fbgn", "log2FoldChange")], by = "fbgn")
cat("Spearman RNA 1 h vs RNA 24 h (all expressed genes):", round(cor(b[, 2], b[, 3], method = "s", use = "complete"), 3), "\n")
