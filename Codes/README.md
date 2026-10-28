# Running the analysis scripts

Run these commands from the repository root. Replace the `input_*` and `output_*` arguments with local files and directories. Run sample-level commands separately for each sample.

## ATAC-seq

```bash
bash Codes/ATAC-seq/process_atac.sh input_reads_1.fq.gz input_reads_2.fq.gz input_dm6_index_prefix output_atac_dir
bash Codes/ATAC-seq/call_peaks_bampe.sh input_sorted.bam output_peak_dir
bash Codes/ATAC-seq/merge_peaks_bedtools.sh input_0h.narrowPeak input_1h.narrowPeak input_3h.narrowPeak input_6h.narrowPeak input_12h.narrowPeak input_24h.narrowPeak output_merged.bed [minimum_score]
Rscript Codes/ATAC-seq/run_atac_wald_1h_3h.R input_counts.csv input_metadata.tsv input_tss_peaks.csv output_wald_dir
python Codes/ATAC-seq/extract_signal.py input_regions.bed input_bigwigs.tsv output_signal.npz
bash Codes/ATAC-seq/select_matching_narrowpeak.sh input_matched.bed input_peaks.narrowPeak output_selected.narrowPeak
```

`minimum_score`, when supplied, filters each narrowPeak input before coordinate merging. The count CSV begins with `peak_id`, followed by sample-count columns. The metadata TSV uses sample IDs as row identifiers and includes `time` and `batch`. The TSS peak CSV includes `peak_id` and may include `matched_gene`. The bigWig manifest TSV contains `sample` and `bigwig` columns.

## RNA-seq

```bash
bash Codes/RNA-seq/download_sra.sh input_accession output_sra_dir
bash Codes/RNA-seq/trim_rna_fastp.sh input_reads.fastq output_trim_dir
bash Codes/RNA-seq/build_star_index.sh input_dm6.fa input_dm6.gtf output_star_index_dir
bash Codes/RNA-seq/align_rna_star.sh input_trimmed.fastq.gz input_star_index_dir input_sample_id output_map_dir
bash Codes/RNA-seq/quantify_rna_salmon.sh input_transcripts.fa input_transcriptome.bam output_quant_dir
Rscript Codes/RNA-seq/run_rna_deseq2.R input_samples.tsv input_tx2gene.tsv output_rna_dir
```

For runs converted as paired files, pass the optional `--split-files` flag to `download_sra.sh`. The RNA sample TSV includes `sample_id`, `path` to `quant.sf`, and `condition` (`ctrl` or `20E`). The transcript map TSV has two headerless columns: transcript ID and gene ID.

## Hi-C, enhancers and GO

```bash
bash Codes/Hi-C/run_juicer_dm6.sh input_juicer_cpu_dir input_dm6.fa input_dm6.chrom.sizes input_dm6_DpnII.txt input_fastq_dir output_project_dir
bash Codes/Hi-C/merge_hic_pc1.sh input_merged30_a.txt input_merged30_b.txt input_dm6.chrom.sizes input_juicer_tools.jar output_hic_dir
Rscript Codes/Hi-C/compartment_qc.R input_positive_counts.txt input_negative_counts.txt input_positive_gc.txt input_negative_gc.txt output_qc.tsv
python Codes/Enhancer/redfly_tss_distance.py input_redfly.tsv input_dm6.gtf output_summary.tsv output_histogram.pdf output_histogram_counts.tsv
Rscript Codes/GO/run_go_quadrants.R input_tss.csv input_enhancer.csv output_go_dir all 3 24
Rscript Codes/GO/run_go_quadrants.R input_tss.csv input_enhancer.csv output_go_dir discordant 1 1
```

Run `run_juicer_dm6.sh` for each Juicer project before combining its `merged30.txt` files. GO input CSVs contain `gene`, `ATAC_testable`, `RNA_testable`, `both_DE`, and `quadrant`; gene IDs are FlyBase IDs. The `all` mode corresponds to 3 h ATAC / 24 h RNA, and the `discordant` mode to 1 h RNA. The hour arguments record the comparison times and label GO plots.

## Software

Command-line tools: Trim Galore, Bowtie2, samtools, MACS2, BEDTools, SRA Toolkit, fastp, STAR, Salmon and Juicer. R packages: DESeq2, tximport, apeglm, ggplot2, clusterProfiler and org.Dm.eg.db. Python packages: pandas, numpy, matplotlib and pyBigWig.

## ATAC–RNA integration

```bash
python Codes/Integration/run_atac_rna_integration.py input_atac.csv input_rna.csv input_gene_mapping.tsv output_integration_dir
```

ATAC CSV columns: `PeakID`, `TSS` (`O` or `X`), `matched_gene`, `enhancer`, `log2FoldChange`, `padj`, and `ClusterLabel`. Gene lists use semicolons. RNA CSV columns: `FBgn`, `log2FoldChange`, and `padj`. Mapping TSV columns: `FLYBASE` and `SYMBOL`. If using the RNA script's output, rename its `gene` column to `FBgn` before supplying it here.

The representative ATAC peak for each gene is the peak with the smallest adjusted p-value. RNA symbols are deduplicated by retaining the first mapped row. Both assays use adjusted p < 0.05 and absolute log2 fold change > log2(1.5). Pearson correlation and linear regression use genes significant in both assays. Direction-specific comparisons use significant RNA genes and matching-direction significant ATAC genes, divided into TSS-only, enhancer-only and shared groups. Pairwise tests use two-sided Mann–Whitney U tests; the supplied method uses unadjusted p-values.

Outputs include representative-peak tables, mapped RNA values, scatter data, correlations, group memberships, group sizes, test results, two scatter PDFs and two violin PDFs. Additional Python dependencies: scipy and adjustText.
