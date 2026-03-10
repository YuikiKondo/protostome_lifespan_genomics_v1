#!/bin/bash
set -euo pipefail

# Input directories (all partial_done dirs)
input_dirs=(
  "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences_partial_done_01"
  "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences_partial_done_02"
  "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences_partial_done_03"
  "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences_partial_done_04"
)

# Output directory
output_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences_split_concatenated"
mkdir -p "$output_dir"

# ------------------------------------------------------------
# 1) Collect all part FASTA files
# ------------------------------------------------------------
tmp_list=$(mktemp)

for d in "${input_dirs[@]}"; do
    if [[ -d "$d" ]]; then
        find "$d" -name "*_part_[0-9][0-9][0-9].fasta"
    else
        echo "Warning: directory not found: $d" >&2
    fi
done > "$tmp_list"

# ------------------------------------------------------------
# 2) Extract unique species prefixes
# ------------------------------------------------------------
sed -E 's#.*/##; s/_part_[0-9]{3}\.fasta$//' "$tmp_list" \
  | sort -u \
  | while read -r prefix; do

    out="${output_dir}/${prefix}_augustus_all_genes_gene_associated_-2kb_+1kb.fasta"

    echo "Merging species: $prefix"

    # --------------------------------------------------------
    # 3) Concatenate parts in numeric order across all dirs
    # --------------------------------------------------------
    grep "/${prefix}_part_" "$tmp_list" \
      | sort \
      | xargs cat > "$out"

    echo "  -> Created: $out"
done

rm -f "$tmp_list"

echo "=== All species merged successfully ==="
