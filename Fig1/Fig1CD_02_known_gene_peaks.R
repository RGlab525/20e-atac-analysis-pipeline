# =============================================================================
# Fig 1C/D (and Supp Fig S6): consensus peaks near 18 known ecdysone-responsive genes and their signal over time
#   Genes : Eip74EF, Eip75B, br (primary early); Hr3, Hr4, Hr39, Eip78C, Blimp-1, Kr-h1; ftz-f1, Eip93F,
#           Eip63E, Eip71CD (MsrA), E23, ImpL2, Cyp18a1, EcR, usp
#   Locus : gene span (org.Dm.eg.db) +/- 10 kb
#   Signal: mean RPGC signal over each consensus peak (Fig1CD_01_peak_signal.sh);
#           change = log2((t + 0.1) / (0 h + 0.1)); background = median and IQR of all 14,003 peaks
#   DAR calls: final DESeq2 table (1 h and 3 h vs 0 h; padj < 0.05, |FC| > 1.5)
# Input : input_consensus_peaks.bed, known_genes/output_all_peaks_mean_signal.tab,
#         input_peak_annotations.xlsx (sheet 2)
# Output: known_genes/output_known_ecdysone_genes_peaks_timecourse.csv, known_genes/output_known_ecdysone_genes_summary.csv,
#         known_genes/output_gene_loci_10kb.tsv
# Source: 261016_revision/08_QC_peak_validity_fragment/scripts/known_genes_02_timecourse.R (J.Y.H., revision R3 Minor 4)
# Changes from the original code: paths via config/config.R; also writes output_gene_loci_10kb.tsv; the check figure is left out
#   (used to choose the genes shown in the Fig 1C tracks; the original file was written separately).
# =============================================================================
suppressMessages({ library(org.Dm.eg.db); library(GenomicRanges); library(readxl) })
options(warn = -1)
source("config/config.R")
proj <- PROJECT_DIR
K    <- out_dir("known_genes")
tp   <- c("0h", "1h", "3h", "6h", "12h", "24h")
genes <- c("Eip74EF", "Eip75B", "br", "Hr3", "Hr4", "Hr39", "Eip78C", "Blimp-1", "Kr-h1",
           "ftz-f1", "Eip93F", "Eip63E", "MsrA", "E23", "ImpL2", "Cyp18a1", "EcR", "usp")
label <- c(MsrA = "Eip71CD (MsrA)")

# gene loci
e  <- mapIds(org.Dm.eg.db, genes, "ENTREZID", "SYMBOL")
a  <- merge(toTable(org.Dm.egCHRLOC), toTable(org.Dm.egCHRLOCEND), by = c("gene_id", "Chromosome"))
a  <- a[a$gene_id %in% e, ]; a$sym <- names(e)[match(a$gene_id, e)]
loc <- do.call(rbind, lapply(split(a, a$sym), function(d) data.frame(sym = d$sym[1], chr = d$Chromosome[1],
          start = min(abs(c(d$start_location, d$end_location))), end = max(abs(c(d$start_location, d$end_location))))))
loc <- loc[match(genes, loc$sym), ]

# peaks + signal + DESeq2
pk  <- read.table(file.path(proj, "input_consensus_peaks.bed"), sep = "\t")
sg  <- read.delim(file.path(K, "output_all_peaks_mean_signal.tab"), check.names = FALSE)
names(sg) <- c("chr", "start", "end", tp)
key <- function(c, s, e) paste(sub("^chr", "", c), s, e)
m   <- as.matrix(sg[match(key(pk$V1, pk$V2, pk$V3), key(sg$chr, sg$start, sg$end)), tp]); m[is.na(m)] <- 0
rownames(m) <- pk$V4
lf  <- log2((m + 0.1) / (m[, "0h"] + 0.1))
ann <- as.data.frame(read_excel(file.path(proj, "input_peak_annotations.xlsx"), sheet = 2))
ann <- ann[match(pk$V4, ann$peak_id), ]
fc  <- log2(1.5)
dar <- function(t, dir) { l <- ann[[paste0("log2FoldChange_", t)]]; p <- ann[[paste0("padj_", t)]]
    !is.na(p) & p < 0.05 & (if (dir > 0) l > fc else l < -fc) }
up1 <- dar("1h", 1); up3 <- dar("3h", 1); dn3 <- dar("3h", -1)
ref <- data.frame(time = factor(tp, levels = tp), med = apply(lf, 2, median),
                  lo = apply(lf, 2, quantile, 0.25), hi = apply(lf, 2, quantile, 0.75))

P  <- GRanges(sub("^chr", "", pk$V1), IRanges(pk$V2 + 1, pk$V3))
hit <- do.call(rbind, lapply(seq_len(nrow(loc)), function(i) {
    w <- GRanges(loc$chr[i], IRanges(loc$start[i] - 10000, loc$end[i] + 10000))
    j <- which(overlapsAny(P, w)); if (!length(j)) return(NULL)
    data.frame(gene = loc$sym[i], idx = j, peak = pk$V4[j], chr = loc$chr[i], start = pk$V2[j], end = pk$V3[j]) }))
# gene loci with >= 1 peak (used to choose the genes shown in the Fig 1C tracks)
lw <- loc[loc$sym %in% hit$gene, ]
write.table(data.frame(gene = lw$sym, chr = lw$chr, start = lw$start - 10000, end = lw$end + 10000),
            file.path(K, "output_gene_loci_10kb.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
hit$class <- ifelse(ann$TSS[hit$idx] == "O", "TSS", "distal")
hit$up_1h <- up1[hit$idx]; hit$up_3h <- up3[hit$idx]; hit$down_3h <- dn3[hit$idx]
hit$log2FC_3h_DESeq2 <- round(ann$log2FoldChange_3h[hit$idx], 2); hit$padj_3h <- signif(ann$padj_3h[hit$idx], 2)
hit <- cbind(hit, round(m[hit$idx, ], 3), setNames(as.data.frame(round(lf[hit$idx, -1], 3)), paste0("log2_", tp[-1], "_vs_0h")))
write.csv(hit[, names(hit) != "idx"], file.path(K, "output_known_ecdysone_genes_peaks_timecourse.csv"), row.names = FALSE)

# per-gene summary
sm <- do.call(rbind, lapply(genes, function(g) { h <- hit[hit$gene == g, ]
    if (!nrow(h)) return(data.frame(gene = g, n_peaks = 0, up_1h = 0, up_3h = 0, down_3h = 0, best_peak = NA,
                                    best_log2FC_3h_DESeq2 = NA, best_padj_3h = NA))
    b <- h[which.max(h$log2FC_3h_DESeq2), ]
    data.frame(gene = g, n_peaks = nrow(h), up_1h = sum(h$up_1h), up_3h = sum(h$up_3h), down_3h = sum(h$down_3h),
               best_peak = b$peak, best_log2FC_3h_DESeq2 = b$log2FC_3h_DESeq2, best_padj_3h = b$padj_3h) }))
sm$gene_label <- ifelse(sm$gene %in% names(label), label[sm$gene], sm$gene)
write.csv(sm, file.path(K, "output_known_ecdysone_genes_summary.csv"), row.names = FALSE)
print(sm, row.names = FALSE)
cat("genes with >=1 peak opening at 3 h (DESeq2):", sum(sm$up_3h > 0), "of", sum(sm$n_peaks > 0), "genes with peaks\n")
