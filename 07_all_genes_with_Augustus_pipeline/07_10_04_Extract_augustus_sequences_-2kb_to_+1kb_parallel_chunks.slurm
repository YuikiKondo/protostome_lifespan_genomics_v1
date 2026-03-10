#!/bin/bash

#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --mem=480G
#SBATCH --time=01-00:00:00
#SBATCH --output=07_06_04_Extract_augustus_sequences_-2kb_to_+1kb_parallel.slurm.log
#SBATCH -p P28

set -euo pipefail

export PATH=$HOME/miniconda3/bin:$PATH
source activate augustus_env

# Define paths
coord_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_coordinates"
split_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_coordinates_split_04"
genome_dir="/home/y-kondo/protostome_lifespan_model/ProtostomeGenomes_unzipped"
partial_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences_partial_04"
partial_done_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences_partial_done_04"
species_list="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/species_in_coordinates_but_missing_in_augustus_sequences.txt"

mkdir -p "$split_dir" "$partial_dir" "$partial_done_dir"

# Step 1: Split coordinate files into 60 equal parts (preserving header)
split_coord_files() {
    sed -n '16,21p' "$species_list" | while read -r species; do
        coord_file="$coord_dir/${species}_augustus_all_genes_gene_associated_-2kb_+1kb.txt"

        if [[ -f "$coord_file" ]]; then
            data_lines=$(tail -n +2 "$coord_file" | wc -l | tr -d ' ')
            if [[ "$data_lines" -gt 0 ]]; then
                echo "Splitting $species with $data_lines lines"

                # If already split, skip
                if ls "$split_dir/${species}_part_"*.txt >/dev/null 2>&1; then
                    echo "Already split: $species (skipping)"
                    continue
                fi

                tmpfile=$(mktemp)
                tail -n +2 "$coord_file" > "$tmpfile"

                split -n l/60 --numeric-suffixes=0 --suffix-length=3 "$tmpfile" "$split_dir/${species}_part_"

                for part_file in "$split_dir/${species}_part_"*; do
                    sed -i '1iRegion\tReference\tStart\tEnd\tStrand\tScore' "$part_file"
                    mv "$part_file" "${part_file}.txt"
                done

                rm -f "$tmpfile"
            else
                echo "Skipping $species: no data lines after header"
            fi
        else
            echo "Warning: Coordinate file not found for $species: $coord_file"
        fi
    done
}

# Step 2: Sequence extraction from split files
process_split_file() {
    local split_file="$1"
    local base species part_suffix fna_file temp_fasta final_fasta

    base=$(basename "$split_file")  # e.g. GCA_..._part_022.txt
    species=$(echo "$base" | sed 's/_part_[0-9]\{3\}\.txt$//')
    part_suffix=$(echo "$base" | grep -o '_part_[0-9]\{3\}')

    # Genome FASTA lookup (try common patterns)
    if [[ -f "${genome_dir}/${species}.fna" ]]; then
        fna_file="${genome_dir}/${species}.fna"
    elif [[ -f "${genome_dir}/${species}" ]]; then
        fna_file="${genome_dir}/${species}"
    else
        echo "Missing genome for species=${species}" >&2
        return 1
    fi

    temp_fasta="${partial_dir}/${species}${part_suffix}.fasta"
    final_fasta="${partial_done_dir}/${species}${part_suffix}.fasta"

    # Skip if already done
    if [[ -s "$final_fasta" ]]; then
        echo "Already done: $final_fasta (skipping)"
        return 0
    fi

    : > "$temp_fasta"

    # Read TSV (skip header) safely
    tail -n +2 "$split_file" | while IFS=$'\t' read -r geneid scaffold start end strand score; do
        seq=$(bedtools getfasta -fi "$fna_file" -bed <(printf "%s\t%s\t%s\n" "$scaffold" "$start" "$end") 2>/dev/null \
              | tail -n +2 | tr -d '\n')

        [[ "$strand" == "-" ]] && seq=$(echo "$seq" | rev | tr 'ACGTacgt' 'TGCAtgca')

        echo ">$geneid augustus $scaffold $start $end $strand $score" >> "$temp_fasta"
        echo "$seq" >> "$temp_fasta"
    done

    # Only move if mv succeeds
    if [[ -s "$temp_fasta" ]]; then
        if mv "$temp_fasta" "$final_fasta"; then
            echo "Finished and moved to $final_fasta"
        else
            echo "ERROR: mv failed for $final_fasta" >&2
            return 1
        fi
    else
        echo "Empty output for $split_file — not moved." >&2
        rm -f "$temp_fasta"
    fi
}

export -f process_split_file
export genome_dir partial_dir partial_done_dir split_dir

# Step 3: Run splitting (only once)
split_coord_files

# Step 4: Run extraction in parallel (only split files from these 5 species)
# (Because split_dir is unique for this run (_04), it will only contain those.)
find "$split_dir" -name "*.txt" | sort \
  | xargs -P 28 -I {} srun --exclusive -n 1 bash -c 'process_split_file "$1"' _ {}
