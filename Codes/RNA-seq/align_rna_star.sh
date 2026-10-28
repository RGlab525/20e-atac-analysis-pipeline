#!/usr/bin/env bash
# STAR single-end two-pass alignment and sorted BAM processing.
set -euo pipefail
if (( $# != 4 )); then
  echo 'Usage: align_rna_star.sh input_trimmed.fastq.gz input_star_index_dir input_sample_id output_dir' >&2
  exit 2
fi
reads=$1; index=$2; sample=$3; outdir=$4
[[ -f "$reads" && -d "$index" ]] || { echo 'FASTQ or STAR index missing' >&2; exit 2; }
mkdir -p "$outdir"
STAR --runThreadN 10 --genomeDir "$index" --readFilesIn "$reads"   --readFilesCommand zcat --outFileNamePrefix "$outdir/output_"   --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts   --outSAMattributes NH HI AS NM MD --outSAMunmapped Within   --outFilterType BySJout --alignSJoverhangMin 8 --outFilterMultimapNmax 20   --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999   --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20   --alignIntronMax 500000 --twopassMode Basic   --outSAMattrRGline ID:"$sample" SM:"$sample" PL:ILLUMINA   > "$outdir/output_star.stdout.log" 2> "$outdir/output_star.stderr.log"
samtools sort -@ 10 -o "$outdir/output_sorted.bam" "$outdir/output_Aligned.out.bam"
samtools index -@ 10 "$outdir/output_sorted.bam"
samtools flagstat -@ 10 "$outdir/output_sorted.bam" > "$outdir/output_sorted.flagstat"
