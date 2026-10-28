# The Dynamic Landscape of Ecdysone-Responsive Chromatin Accessibility

Analysis and existing figure scripts for the 20E_ATAC study, combining code provided by S.C. and J.Y.H.

## Organization

| Folder | Contents |
|---|---|
| `00_preprocessing/` | ATAC trimming/alignment, replicate merging, Tn5 shift, RPGC bigWigs and read summaries |
| `01_peak_processing/` | MACS2 BAMPE peak calling, consensus peak union, featureCounts, count-matrix assembly, representative TSS selection, peak–gene annotation and exact-coordinate peak selection |
| `02_differential_accessibility/` | 1 h / 3 h Wald contrasts and TSS/enhancer fuzzy c-means groups |
| `RNA_seq/` | SRA, fastp, STAR, Salmon, tximport and DESeq2 |
| `HiC/` | Juicer, merged contacts, 250-kb PC1 and compartment QC |
| `Fig1/`–`Fig6/` | Figure-associated analyses and supplied plotting code |
| `SuppFig/` | Sample QC, intergenic/gene-body signal and ATAC–RNA GO analyses |
| `external_data/` | Micro-C, STARR-seq, ChIP-seq, 1 h RNA and small RNA supplementary analyses |
| `config/` | Shared input/output roots, time points and analysis constants |

## Usage

See [RUNNING.md](RUNNING.md) for commands and software, [INPUTS.md](INPUTS.md) for input formats, and [FIGURE_CODE.md](FIGURE_CODE.md) for the Figure-to-code index.

Local data and result names use `input_*` and `output_*` placeholders. Public accessions, reference identifiers, table columns and analysis parameters are retained. Run shared-config scripts from the repository root with `ATAC_PROJECT` and `ATAC_OUT` set to the desired directories.

ATAC coverage uses the supplied RPGC workflow. Figure 2 uses the revised consensus-peak TSS/distal region sets. Figure 6 uses its dedicated peak-annotation input. Existing plots are retained with their analysis scripts; the Figure folders also include scripts that produce statistical tables and track inputs.
