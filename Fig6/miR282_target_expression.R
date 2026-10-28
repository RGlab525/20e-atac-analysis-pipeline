#!/usr/bin/env Rscript
# Plot RNA-seq expression of selected miR-282 target genes.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript miR282_target_expression.R <input_mir282.csv> <output_plot.pdf>")
}
input <- args[[1]]
output <- args[[2]]
if (!file.exists(input)) stop("Input CSV file not found: ", input)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
df <- read.csv(input)
required <- c("con", "Symbol", "log2FC", "padj")
missing <- setdiff(required, names(df))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
df <- df %>% mutate(
  FC_size = abs(log2FC),
  Direction = ifelse(log2FC > 0, "Up", "Down")
)
p <- ggplot(df, aes(x = con, y = Symbol)) +
  geom_tile(aes(fill = padj), color = "white") +
  scale_fill_gradient(low = "grey10", high = "white", trans = "log10", name = "padj") +
  geom_point(aes(size = FC_size, color = Direction)) +
  scale_color_manual(values = c("Up" = "red3", "Down" = "dodgerblue3")) +
  scale_size(range = c(3, 10), name = "|Log2FC|") +
  theme_minimal() +
  labs(title = "miR282 target\ngene expression", x = "", y = "") +
  theme(axis.text = element_text(size = 11, color = "black", face = "bold"),
        panel.grid = element_blank())
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
ggsave(filename = output, width = 2.7, height = 5.5, plot = p)
