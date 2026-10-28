#!/usr/bin/env bash
# =============================================================================
# Analysis code: S.C.; selected source version retained.
# Inputs and outputs: generic input_* / output_* paths supplied as command-line arguments.
# Method and command examples: see RUNNING.md.
# =============================================================================
# Build a dm6 STAR index with sjdbOverhang 49.
set -euo pipefail
if (( $# != 3 )); then
  echo 'Usage: build_star_index.sh input_dm6.fa input_dm6.gtf output_star_index_dir' >&2
  exit 2
fi
[[ -f "$1" && -f "$2" ]] || { echo 'Reference input missing' >&2; exit 2; }
mkdir -p "$3"
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir "$3"   --genomeFastaFiles "$1" --sjdbGTFfile "$2" --sjdbOverhang 49
