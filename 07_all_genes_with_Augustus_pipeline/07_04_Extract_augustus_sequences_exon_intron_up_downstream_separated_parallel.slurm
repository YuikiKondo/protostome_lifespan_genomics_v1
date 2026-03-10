#!/bin/bash

#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --mem=480G
#SBATCH --time=01-00:00:00
#SBATCH --output=/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/07_04_Extract_augustus_sequences_updown_exon_intron_separated_parallel.log
#SBATCH -p P28

export PATH=$HOME/miniconda3/bin:$PATH
source activate augustus_env

# Define paths
coord_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/augustus_coordinates"
genome_dir="/home/y-kondo/protostome_lifespan_model/ProtostomeGenomes_unzipped"
output_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/augustus_sequences"
temp_dir="/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/tmp_augustus_sequences_${SLURM_JOB_ID}_${timestamp}"

mkdir -p "$output_dir" "$temp_dir"
echo "Temporary directory: $temp_dir"

# Function to process a single coordinate file
process_gene_file() {
    coord_file="$1"
    filename=$(basename "$coord_file")
    species=$(echo "$filename" | sed 's/_augustus_all_genes\.txt$//')
    tmp_fasta="${temp_dir}/${species}_augustus_all_genes.fasta"
    final_fasta="${output_dir}/${species}_augustus_all_genes.fasta"
    fna_file="${genome_dir}/${species}.fna"

    echo "Processing $coord_file"

    if [[ ! -f "$fna_file" ]]; then
        echo "Error: Genome file $fna_file not found!"
        return 1
    fi

    > "$tmp_fasta"

    tail -n +2 "$coord_file" | while read -r line; do
        geneid=$(echo "$line" | awk '{print $1}')
        scaffold=$(echo "$line" | awk '{print $2}')
        start=$(echo "$line" | awk '{print $3}')
        end=$(echo "$line" | awk '{print $4}')
        strand=$(echo "$line" | awk '{print $5}')
        score=$(echo "$line" | awk '{print $6}')

        sequence=$(bedtools getfasta -fi "$fna_file" -bed <(echo -e "${scaffold}\t${start}\t${end}") 2>/dev/null | tail -n +2)

        if [[ "$strand" == "-" ]]; then
            sequence=$(echo "$sequence" | rev | tr 'ACGTacgt' 'TGCAtgca')
        fi

        echo ">$geneid augustus $scaffold $start $end $strand $score" >> "$tmp_fasta"
        echo "$sequence" >> "$tmp_fasta"
    done

    # Move the completed file to output_dir
    mv "$tmp_fasta" "$final_fasta"
    echo "Finished writing to $final_fasta"
}

export -f process_gene_file
export coord_dir genome_dir output_dir temp_dir

# Process the 650 coordinate files that haven't been processed yet
find "$coord_dir" -name "*_augustus_all_genes.txt" | sort | sed -n '1,650p' | while read coord_file; do
    species=$(basename "$coord_file" | sed 's/_augustus_all_genes\.txt$//')
    final_fasta="${output_dir}/${species}_augustus_all_genes.fasta"
    if [[ ! -s "$final_fasta" ]]; then
        echo "$coord_file"
    fi
done | xargs -P 28 -I {} srun --exclusive -n 1 bash -c 'process_gene_file "$1"' _ {}

