# =============================================================================
# Supplementary Fig. S2 (intergenic and gene-body signal), step 1: region files
#   Intergenic (S2A): 3-kb regions downstream of transcript ends (input_transcript_downstream_3kb.bed) that do not overlap
#     annotated enhancers, representative TSSs or gene bodies. Original command:
#       bedtools intersect -a input_transcript_downstream_3kb.bed -b input_annotated_enhancers_dm6.bed input_representative_tss.bed input_dm6_transcripts.bed -v \
#         > input_intergenic_reference.bed
#     The same rule is applied here with GenomicRanges (overlap of >= 1 bp excludes a region), and the result
#     is checked against the original file.
#   Gene body (revision, R2-4): transcript from TSS +1 kb to TES (first 1 kb removed so that promoters are not
#     included; remaining length >= 1 kb), excluding annotated enhancers and representative TSSs (same rule as S2A).
#
# Input (relative to ATAC_PROJECT):
#   input_transcript_downstream_3kb.bed           3 kb downstream of each transcript end (UCSC Table Browser, dm6 RefSeq)
#   input_dm6_transcripts.bed           transcript bodies (UCSC Table Browser, dm6 RefSeq)
#   input_annotated_enhancers_dm6.bed   EnhancerAtlas 2.0 D. melanogaster enhancers lifted from dm3 to dm6
#   input_representative_tss.bed representative TSSs
#   input_intergenic_reference.bed   original intergenic regions (only for the check)
#   input_consensus_peaks.bed     consensus peaks (only to report how many gene bodies contain a peak)
# Output ($ATAC_OUT/SuppFig_S2_regions/): output_intergenic_down3000_unique.bed, output_genebody_TSSplus1kb_to_TES_unique.bed
#
# Source: 2025-03-10_ATAC_TSS_enrichment (bedtools command above, line 1267);
#   261016_revision/03_Fig2AB_heatmap_boxplot/scripts/01c_make_genebody_regions.R
# =============================================================================
source("config/config.R")
suppressMessages(library(GenomicRanges)); options(warn = -1)   # hide warnings about contigs missing from one of the files
reg <- out_dir("SuppFig_S2_regions")
rd <- function(f, n = 6) read.table(inp(f), sep = "\t", quote = "", comment.char = "", fill = TRUE)[, 1:n]
gr <- function(x) GRanges(x$V1, IRanges(x$V2 + 1, x$V3))            # BED (0-based) -> 1-based closed
excl <- function(x, files) { k <- rep(FALSE, nrow(x)); for (f in files) k <- k | overlapsAny(gr(x), gr(rd(f, 3))); x[!k, ] }
wr <- function(x, f) write.table(x, file.path(reg, f), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 1) intergenic regions of the original S2A
s2a <- excl(rd("input_transcript_downstream_3kb.bed"), c("input_annotated_enhancers_dm6.bed", "input_representative_tss.bed", "input_dm6_transcripts.bed"))
orig <- rd("input_intergenic_reference.bed")
key <- function(x) paste(x$V1, x$V2, x$V3, x$V4)
cat("intergenic: R", nrow(s2a), "rows / original bedtools", nrow(orig), "rows; identical set:", setequal(key(s2a), key(orig)), "\n")
stopifnot(setequal(key(s2a), key(orig)))
wr(orig, "output_intergenic_down3000_unique.bed")                          # keep the original row order (deepTools input of S2A)

# 2) gene bodies
g <- rd("input_dm6_transcripts.bed")
s <- ifelse(g$V6 == "+", g$V2 + 1000, g$V2); e <- ifelse(g$V6 == "+", g$V3, g$V3 - 1000)
gb <- data.frame(V1 = g$V1, V2 = s, V3 = e, V4 = g$V4, V5 = 0, V6 = g$V6)[e - s >= 1000, ]
cat("transcripts with body >= 1 kb after removing the first 1 kb:", nrow(gb), "\n")
gb <- excl(gb, c("input_annotated_enhancers_dm6.bed", "input_representative_tss.bed"))
gb <- gb[order(gb$V1, gb$V2), ]
cat("after excluding annotated enhancers / representative TSSs:", nrow(gb), "rows,", nrow(unique(gb[, 1:3])), "unique coordinates\n")
pk <- read.table(inp("input_consensus_peaks.bed"), sep = "\t")
u <- unique(gb[, 1:3]); cat("unique gene bodies overlapping a consensus peak:", sum(overlapsAny(gr(u), gr(pk))), "of", nrow(u), "\n")
wr(gb, "output_genebody_TSSplus1kb_to_TES_unique.bed")
