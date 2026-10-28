# =============================================================================
# R1 Major 2 (1 h RNA-seq, GSE304853), step 1: union of all RefSeq exons per gene (dm6), used to turn the authors'
#   RNA-seq bigWig files into gene-level values (a count matrix is not available on GEO)
# Input:  input_dm6_transcripts.bed (UCSC RefSeq BED12, dm6); org.Dm.eg.db (REFSEQ -> FLYBASE, SYMBOL)
# Output: $ATAC_OUT/RNA_1h_GSE304853/output_gene_exons_union_dm6.bed
# Source: 261016_revision/13_R1_Major2_RNA_1h/scripts/01_gene_exons.R (paths only changed)
# =============================================================================
source("config/config.R")
suppressMessages({ library(org.Dm.eg.db); library(GenomicRanges) }); options(warn = -1)
O <- out_dir("RNA_1h_GSE304853")
w <- read.table(inp("input_dm6_transcripts.bed"), sep = "\t", quote = "", comment.char = "")
rs <- toTable(org.Dm.egREFSEQ); fb <- toTable(org.Dm.egFLYBASE); sy <- toTable(org.Dm.egSYMBOL)
eg <- rs$gene_id[match(sub("\\.\\d+$", "", w$V4), rs$accession)]
w$fbgn <- fb$flybase_id[match(eg, fb$gene_id)]; w$sym <- sy$symbol[match(eg, sy$gene_id)]
w <- w[!is.na(w$fbgn), ]
ex <- do.call(rbind, lapply(seq_len(nrow(w)), function(i) {
  sz <- as.integer(strsplit(w$V11[i], ",")[[1]]); st <- as.integer(strsplit(w$V12[i], ",")[[1]])
  data.frame(chr = w$V1[i], start = w$V2[i] + st, end = w$V2[i] + st + sz, fbgn = w$fbgn[i]) }))
g <- GRanges(ex$chr, IRanges(ex$start + 1, ex$end), fbgn = ex$fbgn)
u <- unlist(reduce(split(g, g$fbgn)))
out <- data.frame(chr = as.character(seqnames(u)), start = start(u) - 1, end = end(u), fbgn = names(u))
out$symbol <- w$sym[match(out$fbgn, w$fbgn)]
write.table(out, file.path(O, "output_gene_exons_union_dm6.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("genes", length(unique(out$fbgn)), " exon intervals", nrow(out), "\n")
