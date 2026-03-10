#!/bin/bash
set -euo pipefail

src_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences_split_concatenated"
dst_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences"

mkdir -p "$dst_dir"

shopt -s nullglob

echo "Source dir: $src_dir"
echo "Dest dir  : $dst_dir"
echo "Files in source: $(ls -1 "$src_dir"/*.fasta 2>/dev/null | wc -l)"

moved=0
skipped=0
failed=0

for src in "$src_dir"/*.fasta; do
  base=$(basename "$src")
  dst="$dst_dir/$base"

  if [[ -e "$dst" ]]; then
    echo "SKIP (exists): $base"
    skipped=$((skipped+1))
    continue
  fi

  echo "TRY MOVE: $base"
  if mv "$src" "$dst"; then
    echo "MOVED: $base"
    moved=$((moved+1))
  else
    echo "FAILED: $base" >&2
    failed=$((failed+1))
  fi
done

echo "Done. moved=$moved skipped=$skipped failed=$failed"