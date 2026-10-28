#!/usr/bin/env bash
# Trim Galore paired-end preprocessing, Bowtie2 dm6 alignment and BAM processing.
set -euo pipefail
if (( $# != 4 )); then
  echo 'Usage: process_atac.sh input_reads_1.fq.gz input_reads_2.fq.gz input_bowtie2_index_prefix output_dir' >&2
  exit 2
fi
reads1=$1; reads2=$2; index=$3; outdir=$4
[[ -f "$reads1" && -f "$reads2" ]] || { echo 'Input FASTQ missing' >&2; exit 2; }
mkdir -p "$outdir"
tmpdir=$(mktemp -d "$outdir/.trim.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT
trim_galore --paired --cores 5 --fastqc -o "$tmpdir" "$reads1" "$reads2"
mapfile -t trimmed1 < <(find "$tmpdir" -maxdepth 1 -name '*_val_1.fq.gz')
mapfile -t trimmed2 < <(find "$tmpdir" -maxdepth 1 -name '*_val_2.fq.gz')
(( ${#trimmed1[@]} == 1 && ${#trimmed2[@]} == 1 )) || { echo 'Expected one trimmed FASTQ pair' >&2; exit 1; }
mv "${trimmed1[0]}" "$outdir/output_trimmed_1.fq.gz"
mv "${trimmed2[0]}" "$outdir/output_trimmed_2.fq.gz"
for extension in html zip; do
  mapfile -t reports < <(find "$tmpdir" -maxdepth 1 -name "*fastqc.${extension}" | LC_ALL=C sort)
  for i in "${!reports[@]}"; do
    mv "${reports[$i]}" "$outdir/output_fastqc_$((i+1)).${extension}"
  done
done
mapfile -t trim_reports < <(find "$tmpdir" -maxdepth 1 -name '*trimming_report.txt' | LC_ALL=C sort)
for i in "${!trim_reports[@]}"; do
  mv "${trim_reports[$i]}" "$outdir/output_trimming_report_$((i+1)).txt"
done
bowtie2 -p 7 -x "$index" -1 "$outdir/output_trimmed_1.fq.gz" -2 "$outdir/output_trimmed_2.fq.gz"   | samtools view -F 1804 -q 20 -bS -o "$outdir/output_mapped.bam" -
samtools sort "$outdir/output_mapped.bam" -o "$outdir/output_sorted.bam"
samtools index "$outdir/output_sorted.bam"
