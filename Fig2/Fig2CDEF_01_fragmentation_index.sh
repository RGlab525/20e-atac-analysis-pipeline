#!/bin/bash
# =============================================================================
# Fig 2C-F, step 1: fragmentation index (FI) per 50-kb genomic bin
#   FI = (number of fragments 10-100 bp) / (number of all fragments) in each bin (Methods, "Fragmentation index")
#   Also writes the list of 50-kb bins that contain at least one consensus ATAC peak;
#   only these bins are used in Fig 2C-F.
#
# Input  (relative to ATAC_PROJECT):
#   output_atac_{t}.merge.properly.sorted.bam   merged replicates, Tn5-shifted (alignmentSieve --ATACshift)
#                                               and coordinate-sorted (00_preprocessing); chromosome names without "chr"
#   input_dm6.chrom.sizes
#   input_consensus_peaks.bed               14,003 high-confidence consensus peaks
# Output ($ATAC_OUT/Fig2CDEF_FI/, same layout as the project folder; step 2 reads it with FI_ROOT=$ATAC_OUT/Fig2CDEF_FI):
#   AFI_analysis/output_atac_{t}.fragments.sorted.bed   one line per fragment: chr, start, start+1, fragment length
#   AFI_analysis/output_atac_{t}.apoptotic_counts.bed   fragments 10-100 bp per bin
#   AFI_analysis/output_atac_{t}.total_counts.bed       all fragments per bin
#   AFI_analysis/output_atac_{t}.AFI.bedgraph           FI per bin (0 if the bin has no fragment)
#   dm6_50000_windows.sorted.bed, dm6_50000_windows.sorted.all_time_peaks.intersect.bed
#
# Source: shell history of the analysis (2025-07-05). The file prefix "AFI" (apoptotic fragmentation index)
#   is the working name of the FI. An earlier 10-kb / 180-250-bp version with the same commands was not used.
# Requires: samtools (v1.18), bedtools (v2.31.1), bedGraphToBigWig (optional, for viewing)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
PROJECT_DIR="$(cd "$PROJECT_DIR" && pwd)"
mkdir -p "$OUT_ROOT"
OUT_ROOT="$(cd "$OUT_ROOT" && pwd)"
OUT="$(out_dir Fig2CDEF_FI)"; cd "$OUT"

WINDOW_SIZE=50000
mkdir -p AFI_analysis

# 1) 50-kb windows over the dm6 genome
bedtools makewindows -g "$PROJECT_DIR/input_dm6.chrom.sizes" -w ${WINDOW_SIZE} > output_dm6_${WINDOW_SIZE}_windows_50.bed
sort -k1,1 -k2,2n output_dm6_${WINDOW_SIZE}_windows_50.bed > output_dm6_${WINDOW_SIZE}_windows.sorted.bed

for sample in $TIMES; do
    BAM_FILE="$OUT_ROOT/00_preprocessing/bigwig/output_atac_${sample}.merge.properly.sorted.bam"
    PREFIX="AFI_analysis/output_atac_${sample}"
    echo ">> ${sample}"

    # 2) one record per fragment (read with positive template length), position = leftmost base
    samtools view -F 4 ${BAM_FILE} | \
        awk 'BEGIN{OFS="\t"} $9>0 {print "chr"$3, $4, $4+1, $9}' > ${PREFIX}.fragments.bed
    sort -k1,1 -k2,2n ${PREFIX}.fragments.bed > ${PREFIX}.fragments.sorted.bed

    # 3) short (10-100 bp) and all fragments per window
    awk '$4>=10 && $4<=100' ${PREFIX}.fragments.sorted.bed > ${PREFIX}.apoptotic_fragments.bed
    bedtools map -a output_dm6_${WINDOW_SIZE}_windows.sorted.bed -b ${PREFIX}.apoptotic_fragments.bed -c 1 -o count > ${PREFIX}.apoptotic_counts.bed
    bedtools map -a output_dm6_${WINDOW_SIZE}_windows.sorted.bed -b ${PREFIX}.fragments.sorted.bed   -c 1 -o count > ${PREFIX}.total_counts.bed

    # 4) FI = short / all
    paste <(cut -f 1,2,3 ${PREFIX}.apoptotic_counts.bed) <(cut -f 4 ${PREFIX}.apoptotic_counts.bed) <(cut -f 4 ${PREFIX}.total_counts.bed) | \
        awk 'BEGIN{OFS="\t"} { if ($5 > 0) {afi = $4 / $5} else {afi = 0}; print $1, $2, $3, afi }' > ${PREFIX}.AFI.bedgraph

    # 5) bigWig for genome-browser viewing (not needed for the figure)
    if command -v bedGraphToBigWig > /dev/null; then
        bedGraphToBigWig ${PREFIX}.AFI.bedgraph "$PROJECT_DIR/input_dm6.chrom.sizes" ${PREFIX}.AFI.bw
    fi
done

# 6) windows that overlap at least one consensus peak (used to filter bins in step 2)
bedtools intersect \
  -a output_dm6_${WINDOW_SIZE}_windows.sorted.bed \
  -b "$PROJECT_DIR/input_consensus_peaks.bed" \
  -u > output_dm6_${WINDOW_SIZE}_windows.sorted.all_time_peaks.intersect.bed
wc -l output_dm6_${WINDOW_SIZE}_windows.sorted.all_time_peaks.intersect.bed
