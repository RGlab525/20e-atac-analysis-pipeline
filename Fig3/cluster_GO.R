# =============================================================================
# TSS_cluster_GO_analysis.R
# 목적: cluster별 gene list → GO Term enrichment 분석 + 시각화
# 입력: input_gene_clusters.tsv (col1=gene, col2=cluster, tab-separated)
# 출력: cluster별 CSV / PDF / PNG (4 set)
# 참고: chacha's GO Term pipeline (org.Dm.eg.db, ggplot dotplot style)
# =============================================================================

library(clusterProfiler)
library(org.Dm.eg.db)
library(ggplot2)
library(enrichplot)

# =============================================================================
# ★ 사용자 설정
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) stop("Usage: Rscript cluster_GO.R input_gene_clusters.tsv output_go_dir [TSS|Enhancer]")
INPUT_FILE <- args[[1]]
REGION <- if (length(args) == 3L) args[[3]] else "TSS"
OUTPUT_DIR <- args[[2]]
TOP_N          <- 15       # 상위 GO term 개수
PVALUE_CUTOFF  <- 0.05
PADJ_METHOD    <- "BH"
ONT            <- "BP"

# =============================================================================
# 출력 디렉토리 생성
# =============================================================================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("✅ 결과 저장 위치:", OUTPUT_DIR, "\n\n")

# =============================================================================
# STEP 1: 데이터 로드 및 cluster 분리
# =============================================================================
cat("📂 STEP 1: 데이터 로드 중...\n")

dat_raw <- read.table(INPUT_FILE, header = FALSE, sep = "\t",
                      col.names = c("gene", "cluster"),
                      stringsAsFactors = FALSE,
                      quote = "")          # ← 이 옵션 추가

# gene 이름 정리: "gene1; gene2; gene3" → 세미콜론 기준으로 모두 분리
dat <- do.call(rbind, lapply(seq_len(nrow(dat_raw)), function(i) {
  genes_split <- trimws(strsplit(dat_raw$gene[i], ";")[[1]])
  data.frame(gene    = genes_split,
             cluster = dat_raw$cluster[i],
             stringsAsFactors = FALSE)
}))

# 빈 문자열, NA 제거
dat <- dat[dat$gene != "" & !is.na(dat$gene), ]

cluster_ids <- sort(unique(dat$cluster))
cat("   원본 행 수:", nrow(dat_raw), "\n")
cat("   세미콜론 분리 후 전체 유전자 수:", nrow(dat), "\n")
cat("   클러스터:", paste(cluster_ids, collapse=", "), "\n\n")

for (cid in cluster_ids) {
  cat(paste0(rep("=", 55), collapse=""), "\n")
  cat("🧬 Cluster", cid, "분석 시작\n")
  cat(paste0(rep("=", 55), collapse=""), "\n")

  # ── gene list 추출 + 중복 제거 ─────────────────────────
  genes_raw  <- dat$gene[dat$cluster == cid]
  genes_raw  <- genes_raw[genes_raw != "" & !is.na(genes_raw)]
  genes      <- unique(genes_raw)
  n_removed  <- length(genes_raw) - length(genes)
  cat("   유전자 수:", length(genes),
      if (n_removed > 0) paste0("  (중복 제거: ", n_removed, "개)") else "", "\n")

  # ── GO Enrichment ───────────────────────────────────────
  cat("   GO Enrichment 분석 중...\n")
  ego <- tryCatch(
    enrichGO(
      gene          = genes,
      OrgDb         = org.Dm.eg.db,
      keyType       = "SYMBOL",
      ont           = ONT,
      pAdjustMethod = PADJ_METHOD,
      pvalueCutoff  = PVALUE_CUTOFF,
      readable      = TRUE
    ),
    error = function(e) { cat("   ❌ enrichGO 에러:", conditionMessage(e), "\n"); NULL }
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("   ⚠️  유의한 GO term 없음. 다음 cluster로 넘어갑니다.\n\n")
    next
  }

  ego_df <- as.data.frame(ego)
  cat("   유의한 GO term 수:", nrow(ego_df), "\n")

  # ── 상위 TOP_N 추출 (pvalue 기준) ──────────────────────
  ego_top <- ego_df[order(ego_df$p.adjust), ]
  ego_top <- head(ego_top, TOP_N)

  # Description 줄바꿈 (긴 이름 처리)
  ego_top$Description <- sapply(
    ego_top$Description,
    function(x) paste(strwrap(x, width = 40), collapse = "\n")
  )

  # ── CSV 저장 ────────────────────────────────────────────
  csv_path <- file.path(OUTPUT_DIR, paste0("output_Cluster", cid, "_GOterm_all.csv"))
  write.csv(ego_df, file = csv_path, row.names = FALSE)
  cat("   📄 CSV 저장:", csv_path, "\n")

  # ── ggplot dotplot────────────────────
  p <- ggplot(ego_top,
              aes(x     = FoldEnrichment,
                  y     = reorder(Description, FoldEnrichment),
                  size  = Count,
                  color = p.adjust)) +
    geom_point() +
    labs(
      title = paste0(REGION, " Cluster ", cid, " — GO Term Enrichment (Top ", TOP_N, ")"),
      x     = "Fold Enrichment",
      y     = "GO Term"
    ) +
    scale_color_gradient(low = "indianred1", high = "indianred4") +
    theme_minimal() +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom",
      legend.text     = element_text(size = 7),
      axis.text.y     = element_text(size = 8)
    ) +
    guides(color = guide_colorbar(barwidth = 10, barheight = 1))

  # PDF 저장
  pdf_path <- file.path(OUTPUT_DIR, paste0("output_Cluster", cid, "_GOterm.pdf"))
  ggsave(pdf_path, plot = p, height = 5, width = 5)
  cat("   🖼️  PDF 저장:", pdf_path, "\n")

  # PNG 저장
  png_path <- file.path(OUTPUT_DIR, paste0("output_Cluster", cid, "_GOterm.png"))
  ggsave(png_path, plot = p, height = 5, width = 5, dpi = 300)
  cat("   🖼️  PNG 저장:", png_path, "\n")

  cat("\n")
}

# =============================================================================
# 완료 요약
# =============================================================================
cat(paste0(rep("=", 55), collapse=""), "\n")
cat("✅ 전체 분석 완료!\n")
cat("   결과 위치:", OUTPUT_DIR, "\n")
output_files <- list.files(OUTPUT_DIR, full.names = FALSE)
cat("   생성된 파일:\n")
for (f in output_files) cat("    -", f, "\n")
cat(paste0(rep("=", 55), collapse=""), "\n")
