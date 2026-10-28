# =============================================================================
# STARR-seq enhancers vs distal peaks with increased / decreased accessibility (Supplementary Fig S4D)
#   distal peaks (TSS == "X"); 0 h vs 3 h: increased = padj < 0.05 & log2FC > log2(1.5), decreased = padj < 0.05 &
#   log2FC < -log2(1.5) (Table S7, 1,278 / 634), unchanged = padj >= 0.05 or NA
#   STARR-seq enhancer = summit (lifted to dm6) with P <= 0.001, +/- 250 bp; a peak overlaps if it overlaps an enhancer of
#   either replicate of the condition (untreated, 20E 41 uM 24 h)
#   test: Fisher's exact test (two-sided), increased or decreased vs unchanged peaks; odds ratio with 95% CI
# Input:  input_peak_annotations.xlsx (sheet "260322_1h_3h"),
#         $ATAC_OUT/STARRseq_GSE47691/output_lifted_dm6.rds (STARRseq_GSE47691_01_overlap.R)
# Output: $ATAC_OUT/STARRseq_GSE47691/output_distal_DAR_STARR_odds_ratio.csv
# Source: 261016_revision/19_MicroC_TAD/scripts/05_Fig5FG_alternative_layout.R (odds-ratio part; plotting left out)
# =============================================================================
source("config/config.R")
suppressMessages({ library(GenomicRanges); library(readxl) })
res <- out_dir("STARRseq_GSE47691")

lifted <- readRDS(file.path(res, "output_lifted_dm6.rds"))$starr
win <- lapply(lifted, function(s) { s <- s[s$P <= 1e-3]; resize(s, 501, fix = "center") })
cond <- list(Control = c("noEcd_rep1", "noEcd_rep2"), `20E` = c("ecd_rep1", "ecd_rep2"))

a <- as.data.frame(read_excel(inp("input_peak_annotations.xlsx"), sheet = "260322_1h_3h"))
d <- a[a$TSS == "X", ]; lf <- as.numeric(d$log2FoldChange_3h); pa <- as.numeric(d$padj_3h)
grp <- ifelse(!is.na(pa) & pa < PADJ_CUT & lf > FC_CUT, "Increased", ifelse(!is.na(pa) & pa < PADJ_CUT & lf < -FC_CUT, "Decreased",
       ifelse(is.na(pa) | pa >= PADJ_CUT, "Unchanged", "other")))
print(table(grp))
g <- GRanges(d$chr, IRanges(d$start, d$end))

o <- do.call(rbind, lapply(names(cond), function(cn) {
    h <- Reduce(`|`, lapply(cond[[cn]], function(k) overlapsAny(g, win[[k]])))
    do.call(rbind, lapply(c("Increased", "Decreased"), function(x) {
        k <- grp == x; u <- grp == "Unchanged"
        ft <- fisher.test(matrix(c(sum(h[k]), sum(!h[k]), sum(h[u]), sum(!h[u])), 2))
        data.frame(STARR = cn, group = x, n = sum(k), pct_overlap = 100 * mean(h[k]), pct_overlap_unchanged = 100 * mean(h[u]),
                   OR = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2], P = ft$p.value) })) }))
write.csv(o, file.path(res, "output_distal_DAR_STARR_odds_ratio.csv"), row.names = FALSE)
print(o, digits = 3)
