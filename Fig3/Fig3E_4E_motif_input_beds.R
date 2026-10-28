# Peak BED files for HOMER motif analysis: TSS / enhancer DARs with increased or decreased accessibility (0 h vs 3 h).
#
# Manuscript: Fig 3E (TSS, 131 increased / 1,268 decreased), Fig 4E (enhancer, 1,278 increased / 634 decreased)
# Input : input_peak_annotations.xlsx, sheet 1 (coordinates, 0 h vs 3 h log2FC / padj)
# Output: $ATAC_OUT/motif_beds/{TSS,enhancer}_3h_{up,down}.bed (chr, start, end, peak_id, ., +)
# Source: new script. HOMER's parameter records show the enhancer runs used 260610_enh_3h_{up,down}.bed
#         (run of 2026-06-12); the BED-writing step itself was not kept.
source("config/config.R")
suppressMessages(library(readxl))
out <- out_dir("motif_beds")

tab <- as.data.frame(read_excel(inp("input_peak_annotations.xlsx"), sheet = 1))
tab$log2FoldChange <- as.numeric(tab$log2FoldChange)
tab$padj <- as.numeric(tab$padj)

write_bed <- function(d, name) {
    bed <- data.frame(d$chr, d$start, d$end, d$peak_id, ".", "+")
    bed <- bed[order(bed[[1]], bed[[2]]), ]
    write.table(bed, file.path(out, paste0("output_", name, ".bed")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat(name, nrow(bed), "\n")
}

for (cls in c("TSS", "enhancer")) {
    d <- tab[tab$TSS == if (cls == "TSS") "O" else "X", ]
    write_bed(d[which(d$log2FoldChange >  FC_CUT & d$padj < PADJ_CUT), ], paste0(cls, "_3h_up"))
    write_bed(d[which(d$log2FoldChange < -FC_CUT & d$padj < PADJ_CUT), ], paste0(cls, "_3h_down"))
}
