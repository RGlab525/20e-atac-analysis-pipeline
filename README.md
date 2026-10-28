# The Dynamic Landscape of Ecdysone-Responsive Chromatin Accessibility

## General information

Analysis scripts for the 20E_ATAC study. Script arguments use `input_*` and `output_*` names in place of study-specific file names and paths. [Command examples and input formats](Codes/README.md) are provided separately.

## Data processing

### ATAC-seq

| Script | Description |
|---|---|
| [process_atac.sh](Codes/ATAC-seq/process_atac.sh) | Paired-end Trim Galore, Bowtie2 alignment, BAM filtering and indexing. |
| [call_peaks_bampe.sh](Codes/ATAC-seq/call_peaks_bampe.sh) | MACS2 peak calling in BAMPE mode. |
| [merge_peaks_bedtools.sh](Codes/ATAC-seq/merge_peaks_bedtools.sh) | BEDTools coordinate union across six time points; supports an optional pre-merge narrowPeak score filter. |
| [run_atac_wald_1h_3h.R](Codes/ATAC-seq/run_atac_wald_1h_3h.R) | DESeq2 Wald analysis for 1 h and 3 h against 0 h, with batch in the design. |
| [extract_signal.py](Codes/ATAC-seq/extract_signal.py) | Extract 10 bp binned bigWig signal around BED region centers. |
| [select_matching_narrowpeak.sh](Codes/ATAC-seq/select_matching_narrowpeak.sh) | Retrieve narrowPeak rows matching coordinates in a supplied BED file. |

### RNA-seq

| Script | Description |
|---|---|
| [download_sra.sh](Codes/RNA-seq/download_sra.sh) | Download and convert SRA reads with fasterq-dump. |
| [trim_rna_fastp.sh](Codes/RNA-seq/trim_rna_fastp.sh) | Trim single-end reads with fastp. |
| [build_star_index.sh](Codes/RNA-seq/build_star_index.sh) | Build a dm6 STAR genome index. |
| [align_rna_star.sh](Codes/RNA-seq/align_rna_star.sh) | Two-pass STAR alignment and sorted BAM output. |
| [quantify_rna_salmon.sh](Codes/RNA-seq/quantify_rna_salmon.sh) | Salmon transcript quantification from an alignment BAM. |
| [run_rna_deseq2.R](Codes/RNA-seq/run_rna_deseq2.R) | tximport and DESeq2 differential expression for 20E versus control. |

### Hi-C

| Script | Description |
|---|---|
| [run_juicer_dm6.sh](Codes/Hi-C/run_juicer_dm6.sh) | Run the dm6/DpnII Juicer workflow. |
| [merge_hic_pc1.sh](Codes/Hi-C/merge_hic_pc1.sh) | Combine Juicer contacts, create a `.hic` file and calculate 250 kb PC1 values. |
| [compartment_qc.R](Codes/Hi-C/compartment_qc.R) | Summarize compartment sign and GC-content checks. |

### Enhancer and GO analyses

| Script | Description |
|---|---|
| [redfly_tss_distance.py](Codes/Enhancer/redfly_tss_distance.py) | Measure the distance from REDfly enhancer centers to the nearest same-chromosome transcript TSS. |
| [run_go_quadrants.R](Codes/GO/run_go_quadrants.R) | GO Biological Process enrichment for ATAC/RNA quadrants. |

### ATAC–RNA integration

| Script | Description |
|---|---|
| [run_atac_rna_integration.py](Codes/Integration/run_atac_rna_integration.py) | Map RNA FlyBase IDs to gene symbols, select the lowest-adjusted-p-value ATAC peak per gene, and calculate correlations and direction-specific TSS/enhancer group comparisons. Includes the original scatter and violin outputs. |

TSS assignments use a separately prepared adjusted dm6 TSS BED reference, selecting the TSS of the longest transcript for each gene. The reference is represented as `input_tss_adjusted.bed`; the integration script consumes the resulting gene annotations in the ATAC input table.
