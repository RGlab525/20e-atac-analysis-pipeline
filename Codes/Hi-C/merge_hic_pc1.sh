#!/usr/bin/env bash
# Merge two Juicer merged30 files, build .hic and extract 250 kb PC1.
# Reverse chr3L PC1 after checking its GC-content association.
set -euo pipefail
if (( $# != 5 )); then
  echo 'Usage: merge_hic_pc1.sh input_merged30_a.txt input_merged30_b.txt input_dm6.chrom.sizes input_juicer_tools.jar output_dir' >&2
  exit 2
fi
a=$1; b=$2; sizes=$3; jar=$4; outdir=$5
for file in "$a" "$b" "$sizes" "$jar"; do [[ -f "$file" ]] || { echo "Input missing: $file" >&2; exit 2; }; done
mkdir -p "$outdir"
cat "$a" "$b" | LC_ALL=C sort -k2,2 -k6,6 > "$outdir/output_merged30_sorted.txt"
java -Xmx64g -jar "$jar" pre "$outdir/output_merged30_sorted.txt" "$outdir/output_merged.hic" "$sizes"
for chrom in chr2L chr2R chr3L chr3R chr4 chrX chrY; do
  java -Xmx32g -jar "$jar" eigenvector NONE "$outdir/output_merged.hic" "$chrom" BP 250000 -p "$outdir/output_${chrom}_pc1.txt"
  awk -v chr="$chrom" 'BEGIN{OFS="\t"} {start=(NR-1)*250000; print chr,start,start+250000,$1}'     "$outdir/output_${chrom}_pc1.txt" > "$outdir/output_${chrom}_pc1.bedGraph"
done
awk 'BEGIN{OFS="\t"} {$4=-$4; print}' "$outdir/output_chr3L_pc1.bedGraph" > "$outdir/output_chr3L_pc1_flipped.bedGraph"
cat "$outdir/output_chr2L_pc1.bedGraph" "$outdir/output_chr2R_pc1.bedGraph"     "$outdir/output_chr3L_pc1_flipped.bedGraph" "$outdir/output_chr3R_pc1.bedGraph"     "$outdir/output_chrX_pc1.bedGraph" > "$outdir/output_compartments_250kb.bedGraph"
