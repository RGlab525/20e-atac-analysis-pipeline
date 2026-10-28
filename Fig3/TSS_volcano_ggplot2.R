#!/usr/bin/env Rscript
# =============================================================================
# TSS volcano plot, ggplot2 version (S.C.).
# Input: input_annotated_wald.csv with TSS, log2FoldChange, padj, matched_gene; optional log_p.
#        log_p is the plotted -log10(adjusted P) column from the supplied analysis table.
# Output: output_tss_volcano.pdf/.png and output_tss_volcano_compact.png.
# Method: TSS == O; padj < 0.05 and |log2FC| > 0.585; label 10 lowest-padj DARs.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript TSS_volcano_ggplot2.R input_annotated_wald.csv output_dir")
final_table_with_gene <- read.csv(args[[1]], check.names = FALSE)
required <- c("TSS", "log2FoldChange", "padj", "matched_gene")
if (!all(required %in% names(final_table_with_gene))) stop("Required input columns: ", paste(required, collapse=", "))
if (!"log_p" %in% names(final_table_with_gene)) final_table_with_gene$log_p <- -log10(final_table_with_gene$padj)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(dplyr)
library(ggrepel)

# 1. TSS인 Peak들만 필터링하고 'change' 열 생성
tss_only_table <- final_table_with_gene %>%
  filter(TSS == "O") %>%
  mutate(change = case_when(
    padj < 0.05 & log2FoldChange > 0.585 ~ "Upregulated",
    padj < 0.05 & log2FoldChange < -0.585 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# 2. Y축 최댓값 재계산 (필터링된 데이터 기준)
max_y_tss <- max(tss_only_table$log_p, na.rm = TRUE)

# 3. Volcano Plot 그리기 (tss_only_table 사용)
p <- ggplot(data = tss_only_table, aes(x = log2FoldChange, y = log_p, color = change)) +
  geom_point(alpha = 0.7, size = 2) + # TSS들만 그리므로 점 크기를 살짝 키워도 예쁩니다

  scale_color_manual(name = "Expression Change",
                     values = c("Upregulated" = "#E64B35FF",
                                "Downregulated" = "#0073C2FF",
                                "Not Significant" = "grey")) +

  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +

  # 상위 유전자 이름 표시 (top_genes 데이터도 TSS=="O"인 애들로 미리 뽑아두면 좋습니다)
  geom_text_repel(data = head(tss_only_table %>% filter(change != "Not Significant") %>% arrange(padj), 10),
                  aes(label = matched_gene), # 또는 유전자 이름 열이 있다면 해당 열 사용
                  box.padding = 0.5, max.overlaps = Inf) +

  scale_y_continuous(limits = c(0, max_y_tss + 2), expand = expansion(mult = c(0, 0.05))) +
  theme_light(base_size = 14) +
  labs(title = paste("Volcano Plot: Representative TSS Peaks Only"),
       subtitle = paste("Filtered by TSS == 'O' | Total TSS Peaks:", nrow(tss_only_table)),
       x = expression(paste("Log"[2], " Fold Change")),
       y = expression(paste("-Log"[10], " Adjusted P-value")))

print(p)

ggsave(file.path(output_dir, "output_tss_volcano.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "output_tss_volcano.pdf"), plot = p, width = 10, height = 8)


# 그래프 파일로 저장
ggsave(file.path(output_dir, "output_tss_volcano_compact.png"), width = 8, height = 6)
