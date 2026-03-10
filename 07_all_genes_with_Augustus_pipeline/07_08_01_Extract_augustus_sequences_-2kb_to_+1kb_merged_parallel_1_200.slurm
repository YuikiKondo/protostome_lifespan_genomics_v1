#!/bin/bash
#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=480G
#SBATCH --time=01-00:00:00
#SBATCH --output=/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/07_04_02_02_Extract_augustus_sequences_-2kb_to_+1kb_merged_parallel_1_200.log
#SBATCH -p P64

set -euo pipefail

export PATH=$HOME/miniconda3/bin:$PATH
source activate augustus_env

# INPUT: merged coordinates
coord_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_coordinates_merged"

genome_dir="/home/y-kondo/protostome_lifespan_model/ProtostomeGenomes_unzipped"

# OUTPUT: merged fasta (separate directory is safer)
output_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences_merged"

timestamp=$(date +%Y%m%d_%H%M%S)
temp_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/tmp_augustus_sequences_${SLURM_JOB_ID}_${timestamp}"

mkdir -p "$output_dir" "$temp_dir"
echo "Temporary directory: $temp_dir"

process_merged_coord_file() {
    local coord_file="$1"
    local filename species fna_file tmp_bed tmp_fa final_fa

    filename=$(basename "$coord_file")
    species=$(echo "$filename" | sed 's/_augustus_all_genes_gene_associated_-2kb_+1kb_merged\.txt$//')

    tmp_bed="${temp_dir}/${species}.bed"
    tmp_fa="${temp_dir}/${species}.fasta"
    final_fa="${output_dir}/${species}_augustus_all_genes_gene_associated_-2kb_+1kb_merged.fasta"

    # Genome lookup
    if [[ -f "${genome_dir}/${species}.fna" ]]; then
        fna_file="${genome_dir}/${species}.fna"
    elif [[ -f "${genome_dir}/${species}" ]]; then
        fna_file="${genome_dir}/${species}"
    else
        echo "Error: Genome file not found for species=${species}" >&2
        return 1
    fi

    echo "Processing merged coords: $coord_file"

    # Build BED4 where NAME is pipe-delimited (no spaces) so bedtools parses it safely.
    # Then convert pipes -> spaces to keep your old FASTA header token style (7 tokens).
    #
    # Output header becomes:
    # ><Region> augustus <scaffold> <start> <end> . 0
    #
    tail -n +2 "$coord_file" \
      | awk -F'\t' 'BEGIN{OFS="\t"}{
          region=$1; ref=$2; start=$3; end=$4;
          name=region"|augustus|"ref"|"start"|"end"|.|0";
          print ref, start, end, name;
        }' > "$tmp_bed"

    if [[ ! -s "$tmp_bed" ]]; then
        echo "Warning: no intervals in $coord_file"
        : > "$tmp_fa"
        mv "$tmp_fa" "$final_fa"
        return 0
    fi

    # One bedtools call per species
    bedtools getfasta -fi "$fna_file" -bed "$tmp_bed" -name -tab \
      | awk -F'\t' '{gsub(/\|/," ",$1); print ">"$1"\n"$2}' > "$tmp_fa"

    mv "$tmp_fa" "$final_fa"
    rm -f "$tmp_bed"
    echo "Wrote: $final_fa"
}

export -f process_merged_coord_file
export coord_dir genome_dir output_dir temp_dir

# Process up to 200 files (same pattern you used)
find "$coord_dir" -name "*_augustus_all_genes_gene_associated_-2kb_+1kb_merged.txt" \
  | sort | sed -n '1,200p' \
  | while read -r coord_file; do
      species=$(basename "$coord_file" | sed 's/_augustus_all_genes_gene_associated_-2kb_+1kb_merged\.txt$//')
      final_fa="${output_dir}/${species}_augustus_all_genes_gene_associated_-2kb_+1kb_merged.fasta"
      if [[ ! -s "$final_fa" ]]; then
          echo "$coord_file"
      fi
    done \
  | xargs -P 64 -I {} srun --exclusive -n 1 bash -c 'process_merged_coord_file "$1"' _ {}
