#!/usr/bin/env bash
# Convert SRA reads with fasterq-dump using ten threads.
set -euo pipefail
if (( $# < 2 || $# > 3 )); then
  echo 'Usage: download_sra.sh input_accession output_dir [--split-files]' >&2
  exit 2
fi
accession=$1; outdir=$2
mkdir -p "$outdir"
extra=()
if (( $# == 3 )); then
  [[ "$3" == --split-files ]] || { echo 'Only --split-files is accepted' >&2; exit 2; }
  extra+=(--split-files)
fi
tmpdir=$(mktemp -d "$outdir/.fasterq.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT
fasterq-dump "$accession" "${extra[@]}" --threads 10 --progress -O "$tmpdir"
mapfile -t reads < <(find "$tmpdir" -maxdepth 1 -name '*.fastq' | LC_ALL=C sort)
case ${#reads[@]} in
  1) mv "${reads[0]}" "$outdir/output_reads.fastq" ;;
  2) mv "${reads[0]}" "$outdir/output_reads_1.fastq"
     mv "${reads[1]}" "$outdir/output_reads_2.fastq" ;;
  *) echo 'Expected one or two FASTQ files from fasterq-dump' >&2; exit 1 ;;
esac
