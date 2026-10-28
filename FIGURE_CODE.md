# Figure code index

| Figure / analysis | Code |
|---|---|
| Fig1 | [Fig1CD_01_peak_signal.sh](Fig1/Fig1CD_01_peak_signal.sh) |
| Fig1 | [Fig1CD_02_known_gene_peaks.R](Fig1/Fig1CD_02_known_gene_peaks.R) |
| Fig2 | [Fig2AB_01_region_beds.R](Fig2/Fig2AB_01_region_beds.R) |
| Fig2 | [Fig2AB_02_deepTools_matrices.sh](Fig2/Fig2AB_02_deepTools_matrices.sh) |
| Fig2 | [Fig2CDEF_01_fragmentation_index.sh](Fig2/Fig2CDEF_01_fragmentation_index.sh) |
| Fig2 | [Fig2G_01_boxplot_windows.R](Fig2/Fig2G_01_boxplot_windows.R) |
| Fig2 | [Fig2G_02_signal_stats_SuppFigS2E.R](Fig2/Fig2G_02_signal_stats_SuppFigS2E.R) |
| Fig3 | [Fig3B_S3C_TSS_GSEA.R](Fig3/Fig3B_S3C_TSS_GSEA.R) |
| Fig3 | [Fig3E_4E_HOMER.sh](Fig3/Fig3E_4E_HOMER.sh) |
| Fig3 | [Fig3E_4E_motif_input_beds.R](Fig3/Fig3E_4E_motif_input_beds.R) |
| Fig3 / Fig4 | [cluster_GO.R](Fig3/cluster_GO.R) |
| Fig3 | [TSS_volcano_ggplot2.R](Fig3/TSS_volcano_ggplot2.R) |
| Fig3 | [homer_top_motifs.py](Fig3/homer_top_motifs.py) |
| Fig4 | [Fig4B_Hr4_boxed_peaks.R](Fig4/Fig4B_Hr4_boxed_peaks.R) |
| Fig4 | [REDfly_enhancer_TSS_distance.py](Fig4/REDfly_enhancer_TSS_distance.py) |
| Fig5 | [ATAC_RNA_correlation_and_groups.py](Fig5/ATAC_RNA_correlation_and_groups.py) |
| Fig5 | [Fig5CD_group_genes_vs_RNA_fisher.R](Fig5/Fig5CD_group_genes_vs_RNA_fisher.R) |
| Fig5 | [Fig5F_T1T2_ATAC_RNA_24h_tests.R](Fig5/Fig5F_T1T2_ATAC_RNA_24h_tests.R) |
| Fig6 | [Fig6A_00_mixed_direction_enhancer_peaks.R](Fig6/Fig6A_00_mixed_direction_enhancer_peaks.R) |
| Fig6 | [miR282_target_expression.R](Fig6/miR282_target_expression.R) |
| SuppFig | [ATAC_RNA_GO_quadrants.R](SuppFig/ATAC_RNA_GO_quadrants.R) |
| SuppFig | [SuppFig_S1_sample_QC_metrics.R](SuppFig/SuppFig_S1_sample_QC_metrics.R) |
| SuppFig | [SuppFig_S2_01_intergenic_genebody_regions.R](SuppFig/SuppFig_S2_01_intergenic_genebody_regions.R) |
| SuppFig | [SuppFig_S2_02_intergenic_genebody_fingerprint.sh](SuppFig/SuppFig_S2_02_intergenic_genebody_fingerprint.sh) |
| SuppFig | [SuppFig_S2_03_intergenic_genebody_values.R](SuppFig/SuppFig_S2_03_intergenic_genebody_values.R) |
| external_data | [GSE139316_FLAG_EcR_Usp_at_E1_P_E2.py](external_data/GSE139316_FLAG_EcR_Usp_at_E1_P_E2.py) |
| external_data | [GSE304853_00_download.sh](external_data/GSE304853_00_download.sh) |
| external_data | [GSE304853_01_H3K27ac_EcR_at_peaks.py](external_data/GSE304853_01_H3K27ac_EcR_at_peaks.py) |
| external_data | [GSE304853_02_H3K27ac_EcR_summary.R](external_data/GSE304853_02_H3K27ac_EcR_summary.R) |
| external_data | [MicroC_GSE246517_00_download.sh](external_data/MicroC_GSE246517_00_download.sh) |
| external_data | [MicroC_GSE246517_01_domains.py](external_data/MicroC_GSE246517_01_domains.py) |
| external_data | [MicroC_GSE246517_02_peak_gene_pairs_TADs.py](external_data/MicroC_GSE246517_02_peak_gene_pairs_TADs.py) |
| external_data | [MicroC_GSE246517_03_contacts_OE.py](external_data/MicroC_GSE246517_03_contacts_OE.py) |
| external_data | [RNA1h_GSE304853_01_gene_exons.R](external_data/RNA1h_GSE304853_01_gene_exons.R) |
| external_data | [RNA1h_GSE304853_02_gene_coverage_from_bigwig.py](external_data/RNA1h_GSE304853_02_gene_coverage_from_bigwig.py) |
| external_data | [RNA1h_GSE304853_03_RNA_1h_vs_ATAC.R](external_data/RNA1h_GSE304853_03_RNA_1h_vs_ATAC.R) |
| external_data | [STARRseq_GSE47691_00_download.sh](external_data/STARRseq_GSE47691_00_download.sh) |
| external_data | [STARRseq_GSE47691_01_overlap.R](external_data/STARRseq_GSE47691_01_overlap.R) |
| external_data | [STARRseq_GSE47691_02_distal_DAR_odds_ratio.R](external_data/STARRseq_GSE47691_02_distal_DAR_odds_ratio.R) |
| external_data | [smallRNA_PRJNA632996_00_download.sh](external_data/smallRNA_PRJNA632996_00_download.sh) |
| external_data | [smallRNA_PRJNA632996_01_count_mature_miRNA.py](external_data/smallRNA_PRJNA632996_01_count_mature_miRNA.py) |
| external_data | [smallRNA_PRJNA632996_02_edgeR_miR282.R](external_data/smallRNA_PRJNA632996_02_edgeR_miR282.R) |

TSS / enhancer time-course groups: [TSS c-means](02_differential_accessibility/04_cmeans_TSS_groups.R), [enhancer c-means](02_differential_accessibility/05_cmeans_enhancer_groups.R).

ATAC PCA is retained in the [Wald script](02_differential_accessibility/01_atac_wald_1h_3h.R); RNA dispersion and PCA plots are retained in the [RNA DESeq2 script](RNA_seq/run_rna_deseq2.R).

Individual headers identify J.Y.H. panel assignments, inputs, outputs and analytical methods. REDfly, ATAC–RNA GO and the existing S.C. plots retain their original analysis logic.
