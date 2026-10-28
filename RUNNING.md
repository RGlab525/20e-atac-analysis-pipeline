# Running the analysis scripts

Run these commands from the repository root. Replace the `input_*` and `output_*` arguments with local files and directories. Run sample-level commands separately for each sample.

## ATAC-seq

```bash
bash 00_preprocessing/01_process_atac.sh input_reads_1.fq.gz input_reads_2.fq.gz input_dm6_index_prefix output_atac_dir
bash 01_peak_processing/01_call_peaks_bampe.sh input_sorted.bam output_peak_dir
bash 01_peak_processing/02_merge_peaks_bedtools.sh input_0h.narrowPeak input_1h.narrowPeak input_3h.narrowPeak input_6h.narrowPeak input_12h.narrowPeak input_24h.narrowPeak output_merged.bed [minimum_score]
Rscript 02_differential_accessibility/01_atac_wald_1h_3h.R input_counts.csv input_metadata.tsv input_tss_peaks.csv output_wald_dir
bash 01_peak_processing/03_select_matching_narrowpeak.sh input_matched.bed input_peaks.narrowPeak output_selected.narrowPeak
```

`minimum_score`, when supplied, filters each narrowPeak input before coordinate merging. The count CSV begins with `peak_id`, followed by sample-count columns. The metadata TSV uses sample IDs as row identifiers and includes `time` and `batch`. The TSS peak CSV includes `peak_id` and may include `matched_gene`.

## RNA-seq

```bash
bash RNA_seq/download_sra.sh input_accession output_sra_dir
bash RNA_seq/trim_rna_fastp.sh input_reads.fastq output_trim_dir
bash RNA_seq/build_star_index.sh input_dm6.fa input_dm6.gtf output_star_index_dir
bash RNA_seq/align_rna_star.sh input_trimmed.fastq.gz input_star_index_dir input_sample_id output_map_dir
bash RNA_seq/quantify_rna_salmon.sh input_transcripts.fa input_transcriptome.bam output_quant_dir
Rscript RNA_seq/run_rna_deseq2.R input_samples.tsv input_tx2gene.tsv output_rna_dir
```

For runs converted as paired files, pass the optional `--split-files` flag to `download_sra.sh`. The RNA sample TSV includes `sample_id`, `path` to `quant.sf`, and `condition` (`ctrl` or `20E`). The transcript map TSV has two headerless columns: transcript ID and gene ID.

## Hi-C, enhancers and GO

```bash
bash HiC/run_juicer_dm6.sh input_juicer_cpu_dir input_dm6.fa input_dm6.chrom.sizes input_dm6_DpnII.txt input_fastq_dir output_project_dir
bash HiC/merge_hic_pc1.sh input_merged30_a.txt input_merged30_b.txt input_dm6.chrom.sizes input_juicer_tools.jar output_hic_dir
Rscript HiC/compartment_qc.R input_positive_counts.txt input_negative_counts.txt input_positive_gc.txt input_negative_gc.txt output_qc.tsv
python Fig4/REDfly_enhancer_TSS_distance.py input_redfly.tsv input_dm6.gtf output_summary.tsv output_histogram.pdf output_histogram_counts.tsv
Rscript SuppFig/ATAC_RNA_GO_quadrants.R input_tss.csv input_enhancer.csv output_go_dir all 3 24
Rscript SuppFig/ATAC_RNA_GO_quadrants.R input_tss.csv input_enhancer.csv output_go_dir discordant 1 1
```

Run `run_juicer_dm6.sh` for each Juicer project before combining its `merged30.txt` files. GO input CSVs contain `gene`, `ATAC_testable`, `RNA_testable`, `both_DE`, and `quadrant`; gene IDs are FlyBase IDs. The `all` mode corresponds to 3 h ATAC / 24 h RNA, and the `discordant` mode to 1 h RNA. The hour arguments record the comparison times and label GO plots.

## Software

Command-line tools: Trim Galore, Bowtie2, samtools, MACS2, BEDTools, SRA Toolkit, fastp, STAR, Salmon and Juicer. R packages: DESeq2, tximport, apeglm, ggplot2, clusterProfiler and org.Dm.eg.db. Python packages: pandas, numpy, matplotlib and pyBigWig.

## ATAC–RNA integration

```bash
python Fig5/ATAC_RNA_correlation_and_groups.py input_atac.csv input_rna.csv input_gene_mapping.tsv output_integration_dir
```

ATAC CSV columns: `PeakID`, `TSS` (`O` or `X`), `matched_gene`, `enhancer`, `log2FoldChange`, `padj`, and `ClusterLabel`. Gene lists use semicolons. RNA CSV columns: `FBgn`, `log2FoldChange`, and `padj`. Mapping TSV columns: `FLYBASE` and `SYMBOL`. If using the RNA script's output, rename its `gene` column to `FBgn` before supplying it here.

The representative ATAC peak for each gene is the peak with the smallest adjusted p-value. RNA symbols are deduplicated by retaining the first mapped row. Both assays use adjusted p < 0.05 and absolute log2 fold change > log2(1.5). Pearson correlation and linear regression use genes significant in both assays. Direction-specific comparisons use significant RNA genes and matching-direction significant ATAC genes, divided into TSS-only, enhancer-only and shared groups. Pairwise tests use two-sided Mann–Whitney U tests; the supplied method uses unadjusted p-values.

Outputs include representative-peak tables, mapped RNA values, scatter data, correlations, group memberships, group sizes, test results, two scatter PDFs and two violin PDFs. Additional Python dependencies: scipy and adjustText.

## Shared settings and J.Y.H. Figure analyses

Run from the repository root. `config/config.sh` and `config/config.R` use:

```bash
export ATAC_PROJECT=/path/to/input_data
export ATAC_OUT=/path/to/output_results
export THREADS=8
```

Input files are described in [INPUTS.md](INPUTS.md). Use the same `ATAC_OUT` for successive steps. Generated files retain their `output_*` names when read by later scripts. External-tool standard filenames such as `quant.sf` and `homerResults.html` retain their required names.

### ATAC coverage (RPGC)

Stage the filtered replicate BAMs as `input_bam/input_atac_rep1_{time}.sorted.bam` and `input_bam/input_atac_rep2_{time}.sorted.bam` under `ATAC_PROJECT` (including BAM indexes). Each sample's `01_process_atac.sh` output `output_sorted.bam` is one such filtered BAM.

```bash
bash 00_preprocessing/02_merge_replicates.sh
bash 00_preprocessing/03_shift_and_RPGC_bigwig.sh
bash 00_preprocessing/08_bam_read_summary.sh
Rscript SuppFig/SuppFig_S1_sample_QC_metrics.R
```

The replicate merge is followed by `alignmentSieve --ATACshift`, coordinate sorting and indexing, then `bamCoverage --normalizeUsing RPGC --effectiveGenomeSize 142573017 --ignoreForNormalization chrY chrM --extendReads --ignoreDuplicates`. The default bin size is 50 bp in deepTools 3.5.4. The original normalization exclusion strings are retained exactly. DESeq2 uses count data; small-RNA edgeR retains its TMM normalization and CPM filter.

### Figure analysis order

```bash
# Fig 1
bash Fig1/Fig1CD_01_peak_signal.sh
Rscript Fig1/Fig1CD_02_known_gene_peaks.R

# Fig 2 and Supplementary S2: regions, matrices, per-locus summaries
Rscript Fig2/Fig2AB_01_region_beds.R
Rscript Fig2/Fig2G_01_boxplot_windows.R
bash Fig2/Fig2AB_02_deepTools_matrices.sh
Rscript SuppFig/SuppFig_S2_01_intergenic_genebody_regions.R
bash SuppFig/SuppFig_S2_02_intergenic_genebody_fingerprint.sh
Rscript SuppFig/SuppFig_S2_03_intergenic_genebody_values.R
Rscript Fig2/Fig2G_02_signal_stats_SuppFigS2E.R
bash Fig2/Fig2CDEF_01_fragmentation_index.sh

# Fig 3/4: groups, GSEA and motifs
Rscript 02_differential_accessibility/04_cmeans_TSS_groups.R
Rscript 02_differential_accessibility/05_cmeans_enhancer_groups.R
Rscript Fig3/TSS_volcano_ggplot2.R input_annotated_wald.csv output_volcano_dir
Rscript Fig3/Fig3B_S3C_TSS_GSEA.R
Rscript Fig3/Fig3E_4E_motif_input_beds.R
bash Fig3/Fig3E_4E_HOMER.sh
python Fig3/homer_top_motifs.py output_homer_dir input_group_label output_motifs.tsv
Rscript Fig4/Fig4B_Hr4_boxed_peaks.R

# Fig 5: additional group-level tests
Rscript Fig5/Fig5CD_group_genes_vs_RNA_fisher.R
Rscript Fig5/Fig5F_T1T2_ATAC_RNA_24h_tests.R

# Fig 6
Rscript Fig6/Fig6A_00_mixed_direction_enhancer_peaks.R
Rscript Fig6/miR282_target_expression.R input_mir282.csv output_mir282.pdf
```

`Fig5/ATAC_RNA_correlation_and_groups.py` retains the supplied correlation and violin plots. The revised Fisher and group-comparison analyses use the separate J.Y.H. scripts listed above. `Fig6/miR282_target_expression.R` requires `con`, `Symbol`, `log2FC` and `padj`. The volcano script preserves the selected ggplot2 code and its 0.585 fold-change cutoff; it accepts `log_p` or computes `-log10(padj)` when that plotting column is absent.

### External data / supplementary analyses

```bash
# Micro-C
bash external_data/MicroC_GSE246517_00_download.sh
python external_data/MicroC_GSE246517_01_domains.py
python external_data/MicroC_GSE246517_02_peak_gene_pairs_TADs.py
python external_data/MicroC_GSE246517_03_contacts_OE.py

# STARR-seq and DHS-seq
bash external_data/STARRseq_GSE47691_00_download.sh
Rscript external_data/STARRseq_GSE47691_01_overlap.R
Rscript external_data/STARRseq_GSE47691_02_distal_DAR_odds_ratio.R

# ChIP-seq signal and 1 h RNA coverage
bash external_data/GSE304853_00_download.sh
python external_data/GSE304853_01_H3K27ac_EcR_at_peaks.py
Rscript external_data/GSE304853_02_H3K27ac_EcR_summary.R
python external_data/GSE139316_FLAG_EcR_Usp_at_E1_P_E2.py
Rscript external_data/RNA1h_GSE304853_01_gene_exons.R
python external_data/RNA1h_GSE304853_02_gene_coverage_from_bigwig.py
Rscript external_data/RNA1h_GSE304853_03_RNA_1h_vs_ATAC.R

# Small RNA
bash external_data/smallRNA_PRJNA632996_00_download.sh "$ATAC_PROJECT/input_smallrna"
python external_data/smallRNA_PRJNA632996_01_count_mature_miRNA.py "$ATAC_PROJECT/input_smallrna"
Rscript external_data/smallRNA_PRJNA632996_02_edgeR_miR282.R
```

The external download scripts preserve GEO/SRA accessions and remote filenames; downloaded local files use `input_*` names. ChIP-seq and 1 h RNA scripts analyze the source studies' processed bigWigs. The 1 h RNA response is defined from replicate fold changes, as documented in that script. It is distinct from the 24 h RNA-seq count-based DESeq2 analysis.

Additional software used by J.Y.H.: deepTools 3.5.4, BEDTools 2.31.1, HOMER 5.1, R 4.4 with e1071, GenomicRanges, rtracklayer, Rsamtools, readxl and edgeR; Python with cooler 0.10.4, cooltools 0.7.1 and openpyxl. The selected volcano script also uses dplyr and ggrepel.

## featureCounts and the count matrix

The recorded command is preserved, including `-F SAF -T 14 -p -M`. Supply the six BAMs in time order. A run with merged BAMs produces six merged-time-point columns. For the batch-aware Wald analysis, keep the two replicate/batch count tables separate and assemble their 12 columns:

```bash
bash 01_peak_processing/04_featurecounts.sh input_peaks.saf output_pre_counts.txt input_pre_0h.bam input_pre_1h.bam input_pre_3h.bam input_pre_6h.bam input_pre_12h.bam input_pre_24h.bam
bash 01_peak_processing/04_featurecounts.sh input_peaks.saf output_new_counts.txt input_new_0h.bam input_new_1h.bam input_new_3h.bam input_new_6h.bam input_new_12h.bam input_new_24h.bam
python 01_peak_processing/05_assemble_count_matrix.py output_pre_counts.txt output_new_counts.txt output_matrix_dir
Rscript 02_differential_accessibility/01_atac_wald_1h_3h.R output_matrix_dir/output_counts.csv output_matrix_dir/output_metadata.tsv input_tss_peaks.csv output_wald_dir
```

The SAF is supplied as the original peak annotation; this step does not alter its coordinates. The assembly script copies counts without normalization and writes `pre`/`new` sample metadata. Source output headers record featureCounts 2.0.8 for `pre` and 2.0.1 for `new` and the merged time-point run; retain the corresponding software versions when reproducing those runs.

## Representative TSSs, peak annotation and final workbook

```bash
# 1. Annotate consensus peaks before differential analysis.
python 01_peak_processing/06_annotate_peaks.py input_transcripts.tsv input_consensus_peaks.bed output_annotation_dir

# 2. Use the TSS peak list with the sample-level count matrix.
Rscript 02_differential_accessibility/01_atac_wald_1h_3h.R input_counts.csv input_metadata.tsv output_annotation_dir/output_tss_peaks.csv output_wald_dir

# 3. Join the saved 1 h / 3 h Wald results into the final annotation workbook.
python 01_peak_processing/06_annotate_peaks.py input_transcripts.tsv input_consensus_peaks.bed output_annotation_dir --wald-1h output_wald_dir/output_wald_1h_all.csv --wald-3h output_wald_dir/output_wald_3h_all.csv
```

The representative TSS is transcript start minus one on the plus strand, transcript end minus one on the minus strand. TSS/distal classification uses 1-bp TSS BED overlaps. Gene assignment to TSS-class peaks retains the original R inclusive matching on the stored coordinate numbers; it is kept separate from BED classification. Distal peaks link to all active representative TSSs in a 10-kb interval window. Active TSSs are those overlapping a consensus peak.

Use `output_tss_representative.bed` as `input_representative_tss.bed` and `output_peak_annotations.xlsx` as `input_peak_annotations.xlsx` for the corresponding downstream Figure scripts.

## TSS / enhancer cluster GO

```bash
Rscript Fig3/cluster_GO.R input_tss_gene_clusters.tsv output_tss_go_dir TSS
Rscript Fig3/cluster_GO.R input_enhancer_gene_clusters.tsv output_enhancer_go_dir Enhancer
```

The supplied method uses `enrichGO`, `org.Dm.eg.db`, `keyType="SYMBOL"`, `ont="BP"`, BH adjustment, `pvalueCutoff=0.05`, and `readable=TRUE`. It retains the source defaults for the background, gene-set size limits and q-value cutoff. It plots the first 15 returned terms ordered by adjusted P, with FoldEnrichment, Count and adjusted P. This is separate from TSS GSEA and the ATAC–RNA quadrant GO analysis. Groups with no significant terms do not produce a plot.

The supplied cluster GO results were checked with clusterProfiler 4.18.4 and org.Dm.eg.db 3.22.0.
