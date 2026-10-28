# Consensus peaks inside the dashed boxes of the Hr4 locus (Fig 4B) with their 1 h / 3 h DESeq2 results.
#
# Manuscript: Fig 4B and its legend (TSS box: 1 TSS peak; enhancer box: 7 separate enhancer peaks),
#             response to Reviewer 1 Minor 3
# Boxes : drawn in Illustrator on the UCSC view chrX:1,931,247-2,010,077; converted to genome coordinates
#         from their position on the figure scale bar: TSS box chrX:1,941,568-1,943,064,
#         enhancer box chrX:1,943,962-1,953,535.
# Input : input_consensus_peaks.bed, input_peak_annotations.xlsx (sheet 2)
# Output: $ATAC_OUT/Fig4B/output_Fig4B_Hr4_boxed_peaks.csv
# Source: 261016_revision/05_Fig4B_Hr4_classification/scripts/hr4_and_classification.R (part 17)
# The track panel itself is a UCSC Genome Browser view.
source("config/config.R")
suppressMessages({ library(GenomicRanges); library(readxl) })
out <- out_dir("Fig4B")

pk  <- read.table(inp("input_consensus_peaks.bed"), sep = "\t")
pkg <- GRanges(pk$V1, IRanges(pk$V2 + 1, pk$V3), name = pk$V4)
ann <- as.data.frame(read_excel(inp("input_peak_annotations.xlsx"), sheet = 2))
a   <- ann[match(pkg$name, ann$peak_id), ]
pkg$class <- ifelse(a$TSS == "O", "TSS", "Enhancer")

boxes <- data.frame(box = c("TSS box", "Enhancer box"), start = c(1941568, 1943962), end = c(1943064, 1953535))
h <- findOverlaps(pkg, GRanges("chrX", IRanges(boxes$start, boxes$end)))
q <- queryHits(h)
hr <- data.frame(box = boxes$box[subjectHits(h)], peak = pkg$name[q], start = start(pkg)[q], end = end(pkg)[q],
                 width_bp = width(pkg)[q], class = pkg$class[q],
                 log2FC_1h = round(a$log2FoldChange_1h[q], 2), padj_1h = signif(a$padj_1h[q], 2),
                 log2FC_3h = round(a$log2FoldChange_3h[q], 2), padj_3h = signif(a$padj_3h[q], 2),
                 linked_gene = a$enhancer[q], tss_gene = a$matched_gene[q])
hr$gap_to_next_peak_bp <- c(diff(hr$start) - head(hr$width_bp, -1), NA)
print(hr, row.names = FALSE)
write.csv(hr, file.path(out, "output_Fig4B_Hr4_boxed_peaks.csv"), row.names = FALSE)
