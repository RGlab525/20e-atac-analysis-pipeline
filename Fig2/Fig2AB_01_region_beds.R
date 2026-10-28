# =============================================================================
# Fig 2A/B, step 1: region files for the TSS and putative regulatory region (distal peak) heatmaps
#   output_A_TSS_peaks_4774.bed     one representative TSS (with strand) for each of the 4,774 TSS-associated consensus peaks;
#                            representative TSS = TSS in input_representative_tss.bed that lies inside the peak,
#                            preferring the TSS of the gene assigned to the peak (matched_gene)
#   output_B_distal_peaks_9229.bed  the 9,229 consensus peaks that do not overlap a representative TSS
#                            (used with --referencePoint center)
#
# Input  (relative to ATAC_PROJECT):
#   input_consensus_peaks.bed                          14,003 consensus peaks
#   input_peak_annotations.xlsx   sheet 2: peak_id, TSS (O/X), matched_gene
#   input_representative_tss.bed                            representative TSS per gene (longest isoform)
# Output ($ATAC_OUT/Fig2AB_regions/): output_A_TSS_peaks_4774.bed, output_B_distal_peaks_9229.bed
#
# Source: 261016_revision/03_Fig2AB_heatmap_boxplot/scripts/make_beds.R (revision, R1 Major 1-B / R3 Minor 4).
#   The submitted Fig 2A/B used all UCSC TSSs (dm6_ucsc_tss.bed) and EnhancerAtlas enhancers
#   (input_annotated_enhancers_dm6.bed); the revised panels use the consensus peaks.
# =============================================================================
source("config/config.R")
suppressMessages({ library(GenomicRanges); library(readxl) })
out <- out_dir("Fig2AB_regions")

pk  <- read.table(inp("input_consensus_peaks.bed"), sep = "\t")
pkg <- GRanges(pk$V1, IRanges(pk$V2 + 1, pk$V3), name = pk$V4)
ann <- as.data.frame(read_excel(inp("input_peak_annotations.xlsx"), sheet = 2))
a   <- ann[match(pkg$name, ann$peak_id), ]
isT <- a$TSS == "O"

x   <- read.table(inp("input_representative_tss.bed"), sep = "\t", quote = "", comment.char = "")
tss <- GRanges(x$V1, IRanges(x$V2 + 1, width = 1), strand = x$V6, gene = x$V4)

tp  <- pkg[isT]; mg <- a$matched_gene[isT]
h   <- suppressWarnings(findOverlaps(tp, tss, ignore.strand = TRUE))
cand <- data.frame(q = queryHits(h), s = subjectHits(h), gene = tss$gene[subjectHits(h)])
cand$match <- mapply(function(g, m) g %in% trimws(strsplit(as.character(m), "[;,]")[[1]]), cand$gene, mg[cand$q])
cand <- cand[order(cand$q, !cand$match), ]
pick <- cand[!duplicated(cand$q), ]
cat(sprintf("TSS peaks %d; with representative TSS %d; matched to peak gene %d\n", length(tp), nrow(pick), sum(pick$match)))
sel <- tss[pick$s]
A <- data.frame(chr = as.character(seqnames(sel)), start = start(sel) - 1, end = start(sel),
                name = paste(tp$name[pick$q], sel$gene, sep = "|"), score = 0, strand = as.character(strand(sel)))
A <- A[order(A$chr, A$start), ]
write.table(A, file.path(out, "output_A_TSS_peaks_4774.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

dp <- pkg[!isT]
B <- data.frame(chr = as.character(seqnames(dp)), start = start(dp) - 1, end = end(dp), name = dp$name, score = 0, strand = ".")
B <- B[order(B$chr, B$start), ]
write.table(B, file.path(out, "output_B_distal_peaks_9229.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("A rows %d, B rows %d\n", nrow(A), nrow(B)))
