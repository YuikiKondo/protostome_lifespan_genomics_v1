#!/bin/bash
#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/07_04.01_calc_total_length_of_augustus_sequence_-2kb_to_+1kb_parallel.slurm.log
#SBATCH -p P64

set -euo pipefail

fasta_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences_merged"
out_csv="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/gene_associated_-2kb_+1kb_merged_total_bp_per_species.csv"

# Write header
printf "Species,Total_bp_gene_associated_region_including_N,Num_merged_intervals\n" > "$out_csv"

# Fail early if no files match
n_files=$(find "$fasta_dir" -maxdepth 1 -type f -name "*_genomic_augustus_all_genes_gene_associated_-2kb_+1kb_merged.fasta" | wc -l)
if [[ "$n_files" -eq 0 ]]; then
  echo "ERROR: No FASTA files found in: $fasta_dir" >&2
  echo "Check filename pattern or directory." >&2
  exit 1
fi

# Parallel processing
find "$fasta_dir" -maxdepth 1 -type f \
  -name "*_genomic_augustus_all_genes_gene_associated_-2kb_+1kb_merged.fasta" \
  | sort \
  | xargs -P 64 -I {} bash -c '
      set -euo pipefail
      fasta="$1"

      species=$(basename "$fasta" | sed "s/_genomic_augustus_all_genes_gene_associated_-2kb_+1kb_merged\.fasta$//")

      # Total bp: sum lengths of all non-header lines
      total_bp=$(awk "
        BEGIN{sum=0}
        /^>/{next}
        {gsub(/[ \t\r\n]/, \"\"); sum += length(\$0)}
        END{print sum}
      " "$fasta")

      # Number of merged intervals: count headers
      nrec=$(awk "/^>/{c++} END{print c+0}" "$fasta")

      # Output ONE clean CSV line (no extra quotes)
      printf "%s,%s,%s\n" "$species" "$total_bp" "$nrec"
  ' _ {} >> "$out_csv"

echo "Finished. Output written to: $out_csv"
