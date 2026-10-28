#!/usr/bin/env Rscript
# GO Biological Process enrichment for ATAC/RNA quadrants.
# The all mode uses four quadrants; discordant uses +- and -+.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6L) {
  stop("Usage: Rscript run_go_quadrants.R input_tss.csv input_enhancer.csv output_dir <all|discordant> <ATAC hour> <RNA hour>")
}
input_paths <- c(TSS=args[[1]], Enhancer=args[[2]])
output_dir <- args[[3]]
mode <- args[[4]]
atac_hour <- args[[5]]
rna_hour <- args[[6]]
if (!mode %in% c("all", "discordant")) stop("Mode must be all or discordant")
for (path in input_paths) if (!file.exists(path)) stop("Input missing: ", path)
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Dm.eg.db)
  library(ggplot2)
})
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
quadrants <- if (mode == "all") c("++", "+-", "-+", "--") else c("+-", "-+")
q_code <- c("++"="pp", "+-"="pn", "-+"="np", "--"="nn")
fdr_cutoff <- 0.05
top_n <- 15L

ratio_to_numeric <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed=TRUE), function(parts) {
    if (length(parts) != 2L) return(NA_real_)
    as.numeric(parts[1]) / as.numeric(parts[2])
  }, numeric(1))
}

save_plot <- function(top, title, prefix, input_n) {
  if (nrow(top)) {
    top$Description <- vapply(top$Description, function(x)
      paste(strwrap(x, width=33), collapse="\n"), character(1))
    plot <- ggplot(top, aes(x=FoldEnrichment,
                            y=reorder(Description, FoldEnrichment),
                            size=Count, color=pvalue)) +
      geom_point() +
      labs(title=title, x="Fold Enrichment", y="GO Term") +
      scale_color_gradient(low="indianred1", high="indianred4") +
      theme_minimal() +
      theme(plot.title=element_text(hjust=0.5, size=12, face="bold"),
            legend.position="bottom", axis.text.y=element_text(size=8))
  } else {
    label <- paste0("No GO BP terms pass BH-adjusted p < 0.05",
                    if (mode == "discordant") " and Count >= 2" else "",
                    "\n(input genes: ", input_n, ")")
    plot <- ggplot() + annotate("text", x=0, y=0, label=label, size=4) +
      xlim(-1,1) + ylim(-1,1) + labs(title=title) + theme_void()
  }
  height <- if (mode == "discordant" && atac_hour == "1") 6.5 else 5
  ggsave(paste0(prefix, ".pdf"), plot=plot, width=5, height=height)
  ggsave(paste0(prefix, ".png"), plot=plot, width=5, height=height, dpi=300)
}

summary_rows <- list()
for (region in names(input_paths)) {
  dat <- read.csv(input_paths[[region]], check.names=FALSE,
                  stringsAsFactors=FALSE)
  required <- c("gene", "ATAC_testable", "RNA_testable", "both_DE", "quadrant")
  if (!all(required %in% names(dat))) stop("Missing gene-level columns in ", input_paths[[region]])
  for (flag in c("ATAC_testable", "RNA_testable", "both_DE"))
    dat[[flag]] <- toupper(as.character(dat[[flag]])) == "TRUE"
  background <- unique(dat$gene[dat$ATAC_testable & dat$RNA_testable])
  background <- background[!is.na(background) & nzchar(background)]
  region_prefix <- paste0("output_", tolower(region))
  write.csv(data.frame(gene=background),
            file.path(output_dir, paste0(region_prefix, "_background.csv")),
            row.names=FALSE)
  for (q in quadrants) {
    genes <- unique(dat$gene[dat$both_DE & dat$quadrant == q])
    genes <- intersect(genes[!is.na(genes)], background)
    prefix <- file.path(output_dir, paste0(region_prefix, "_", q_code[[q]]))
    write.csv(data.frame(gene=genes), paste0(prefix, "_input_genes.csv"),
              row.names=FALSE)
    ego <- if (length(genes)) suppressMessages(enrichGO(
      gene=genes, universe=background, OrgDb=org.Dm.eg.db,
      keyType="FLYBASE", ont="BP", pAdjustMethod="BH",
      pvalueCutoff=1, qvalueCutoff=1, minGSSize=10,
      maxGSSize=500, readable=TRUE)) else NULL
    raw <- if (is.null(ego)) data.frame() else as.data.frame(ego)
    if (nrow(raw) && !"FoldEnrichment" %in% names(raw)) {
      raw$FoldEnrichment <- ratio_to_numeric(raw$GeneRatio) /
                            ratio_to_numeric(raw$BgRatio)
    }
    write.csv(raw, paste0(prefix, "_go_all.csv"), row.names=FALSE)
    significant <- if (!nrow(raw)) raw else raw[
      !is.na(raw$p.adjust) & raw$p.adjust < fdr_cutoff &
        is.finite(raw$FoldEnrichment) &
        (mode == "all" | raw$Count >= 2), , drop=FALSE]
    if (nrow(significant)) {
      significant <- significant[order(significant$p.adjust,
                                       significant$pvalue,
                                       -significant$FoldEnrichment), , drop=FALSE]
    }
    write.csv(significant, paste0(prefix, "_go_significant.csv"),
              row.names=FALSE)
    top <- head(significant, top_n)
    title <- paste0("ATAC ", atac_hour, " h / RNA ", rna_hour, " h: ",
                    region, " ", q, "\nGO BP Top ", top_n)
    save_plot(top, title, paste0(prefix, "_go_plot"), length(genes))
    summary_rows[[length(summary_rows)+1L]] <- data.frame(
      region=region, quadrant=q, input_genes=length(genes),
      background_genes=length(background), GO_terms_tested=nrow(raw),
      GO_terms_FDR_0_05=nrow(significant), GO_terms_plotted=nrow(top))
  }
}
summary <- do.call(rbind, summary_rows)
write.csv(summary, file.path(output_dir, "output_go_summary.csv"),
          row.names=FALSE)
writeLines(c(
  paste("ATAC hour:", atac_hour, "RNA hour:", rna_hour, "mode:", mode),
  "Universe: region-specific genes testable in both assays.",
  "GO BP, FLYBASE IDs, BH correction, term size 10-500, FDR < 0.05.",
  if (mode == "discordant") "Discordant mode also requires Count >= 2." else
    "All mode includes all four ATAC/RNA direction quadrants.",
  paste("org.Dm.eg.db:", as.character(packageVersion("org.Dm.eg.db"))),
  paste("clusterProfiler:", as.character(packageVersion("clusterProfiler")))
), file.path(output_dir, "output_go_methods.txt"))
print(summary)
