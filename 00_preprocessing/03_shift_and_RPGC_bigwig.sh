#!/usr/bin/env bash
# =============================================================================
# Tn5 shift and RPGC-normalized coverage tracks (one bigWig per time point, merged replicates)
#   Input : $ATAC_OUT/00_preprocessing/bam/output_atac_{t}.merge.sorted.bam
#   Output: output_atac_{t}.merge.properly.sorted.bam  (shifted, sorted; also used for the FI analysis, Fig2)
#           output_atac_{t}.tss.bw                     (RPGC, 50-bp bins = bamCoverage default)
#   The bigWig files are used for all genome-browser views (Fig 1C, 4B, 6B, Supp S3/S4),
#   deepTools heatmaps/profiles (Fig 2A/B) and peak-signal summaries.
# Source: 2025-03-10_ATAC_TSS_enrichment (command log, J.Y.H.); parameters copied unchanged.
# Note: chromosome names in the BAMs are "Y" and "M", so "--ignoreForNormalization chrY chrM" may not
#       have matched any chromosome; the option is kept as it was run.
# Tools: deepTools v3.5.4 (alignmentSieve, bamCoverage), samtools
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT=$(out_dir 00_preprocessing/bigwig)

for t in $TIMES; do
    alignmentSieve \
      -b "$OUT_ROOT/00_preprocessing/bam/output_atac_${t}.merge.sorted.bam" \
      -o "$OUT/output_atac_${t}.merge.sorted.shift.bam" \
      --ATACshift \
      --verbose

    samtools sort "$OUT/output_atac_${t}.merge.sorted.shift.bam" -o "$OUT/output_atac_${t}.merge.properly.sorted.bam"
    samtools index "$OUT/output_atac_${t}.merge.properly.sorted.bam"

    bamCoverage \
      -b "$OUT/output_atac_${t}.merge.properly.sorted.bam" \
      -o "$OUT/output_atac_${t}.tss.bw" \
      --normalizeUsing RPGC \
      --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" \
      --ignoreForNormalization chrY chrM \
      --extendReads \
      --ignoreDuplicates \
      -p "$THREADS"
done
