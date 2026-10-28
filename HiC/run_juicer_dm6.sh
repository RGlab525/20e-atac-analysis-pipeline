#!/usr/bin/env bash
# =============================================================================
# Analysis code: S.C.; selected source version retained.
# Inputs and outputs: generic input_* / output_* paths supplied as command-line arguments.
# Method and command examples: see RUNNING.md.
# =============================================================================
# Run once per Juicer project to create its merged30.txt contact file.
set -euo pipefail
if (( $# != 6 )); then
  echo 'Usage: run_juicer_dm6.sh input_juicer_cpu_dir input_dm6.fa input_dm6.chrom.sizes input_dm6_DpnII.txt input_fastq_dir output_project_dir' >&2
  exit 2
fi
cpu=$1; fasta=$2; sizes=$3; sites=$4; fastq=$5; project=$6
for file in "$fasta" "$sizes" "$sites"; do [[ -f "$file" ]] || { echo "Input missing: $file" >&2; exit 2; }; done
[[ -x "$cpu/juicer.sh" && -d "$fastq" ]] || { echo 'Juicer executable or FASTQ directory missing' >&2; exit 2; }
mkdir -p "$project"
if [[ -e "$project/fastq" || -L "$project/fastq" ]]; then
  [[ "$(realpath "$project/fastq")" == "$(realpath "$fastq")" ]] || {
    echo 'Existing project/fastq points to a different input directory' >&2
    exit 2
  }
else
  ln -s "$(realpath "$fastq")" "$project/fastq"
fi
"$cpu/juicer.sh" -g dm6 -s DpnII -d "$project" -D "$cpu" -z "$fasta" -p "$sizes" -y "$sites"
