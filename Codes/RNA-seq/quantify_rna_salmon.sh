#!/usr/bin/env bash
# Alignment-based Salmon quantification with SR library type and bias correction.
set -euo pipefail
if (( $# != 3 )); then
  echo 'Usage: quantify_rna_salmon.sh input_transcripts.fa input_transcriptome.bam output_quant_dir' >&2
  exit 2
fi
[[ -f "$1" && -f "$2" ]] || { echo 'Transcript FASTA or BAM missing' >&2; exit 2; }
salmon quant -t "$1" --libType SR -a "$2" -o "$3"   --fldMean 250 --fldSD 80 --gcBias --seqBias
# quant.sf is written inside output_quant_dir.
