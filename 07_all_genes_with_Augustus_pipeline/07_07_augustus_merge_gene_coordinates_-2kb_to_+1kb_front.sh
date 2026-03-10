#!/bin/bash
set -euo pipefail

export PATH=$HOME/miniconda3/bin:$PATH
source activate augustus_env   # bedtools is available here

in_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_coordinates"
out_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_coordinates_merged"
mkdir -p "$out_dir"

merge_one() {
  local coord_file="$1"
  local base species out tmp

  base=$(basename "$coord_file")
  species=$(echo "$base" | sed 's/_augustus_all_genes_gene_associated_-2kb_+1kb\.txt$//')
  out="${out_dir}/${species}_augustus_all_genes_gene_associated_-2kb_+1kb_merged.txt"
  tmp=$(mktemp)

  # Convert TSV -> BED3 (Reference, Start, End), sort, merge overlaps/touching
  tail -n +2 "$coord_file" \
    | awk -F'\t' 'BEGIN{OFS="\t"}{
        ref=$2; s=$3; e=$4;
        if (s < 0) s=0;
        if (e > s) print ref, s, e;
      }' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > "$tmp"

  # Write merged TSV with header
  {
    echo -e "Region\tReference\tStart\tEnd\tStrand\tScore"
    awk 'BEGIN{OFS="\t"}{
        print "gene_associated_-2kb_+1kb_merged_"NR, $1, $2, $3, ".", 0
      }' "$tmp"
  } > "$out"

  rm -f "$tmp"
  echo "Merged: $coord_file -> $out"
}

export -f merge_one
export out_dir

find "$in_dir" -name "*_augustus_all_genes_gene_associated_-2kb_+1kb.txt" \
  | sort \
  | xargs -P 28 -I {} bash -c 'merge_one "$1"' _ {}
