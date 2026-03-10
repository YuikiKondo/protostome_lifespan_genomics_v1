#!/bin/bash
## It is recommended that around 200 files are processed in one job.

#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1                       # Request nodes
#SBATCH --ntasks=64                     # Total tasks across nodes
#SBATCH --cpus-per-task=1               # One CPU per task
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00-24:00:00
#SBATCH --output=/home/username/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/05_03_Extract_busco_sequences_updown_exon_intron_separted_parallel.log
#SBATCH -p P128

# Load necessary modules or export paths (depending on your environment setup)
export PATH=$HOME/miniconda3/bin:$PATH

# Activate the conda environment with BUSCO
source activate busco_env

# Define paths
busco_dir="/home/username/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/busco_coordinates_including_dup/"
genome_dir="/home/username/protostome_lifespan_model/ProtostomeGenomes_unzipped"   # Directory with genome files
output_dir="/home/username/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/busco_sequence1_650/"   # Directory for output FASTA files

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Define a function to process a single gene file
process_gene_file() {
    gene_file="$1"
    gene_filename=$(basename "$gene_file")

    echo "Processing gene file: $gene_file"

    # Extract the species name from the gene file name
    species=$(basename "$gene_file" | sed -E 's/^BUSCO_//' | sed -E 's/_(complete|fragmented|duplicated)_genes\.txt$//')

    # Determine gene type from the filename
    if [[ "$gene_file" =~ (complete|fragmented|duplicated) ]]; then
        gene_type="${BASH_REMATCH[1]}"
    else
        echo "Error: Gene type not found in file name"
        return
    fi

    # Output file for this species and gene type
    output_fasta="${output_dir}/${species}_${gene_type}_genes.fasta"

    # Get the corresponding .fna file for the species
    fna_file="${genome_dir}/${species}"
    
    # Skip if the .fna file doesn't exist
    if [[ ! -f "$fna_file" ]]; then
        echo "Error: Genomic file for species $species does not exist!"
        return
    fi

    # Clear or create the output file for this gene type
    > "$output_fasta"

    # Process each line in the gene file, skipping the header by reading from the second line
    tail -n +2 "$gene_file" | while read -r line; do
        geneid=$(echo "$line" | awk '{print $1}')
        scaffold=$(echo "$line" | awk '{print $2}')
        start=$(echo "$line" | awk '{print $3}')
        end=$(echo "$line" | awk '{print $4}')
        strand=$(echo "$line" | awk '{print $5}')
        score=$(echo "$line" | awk '{print $6}')

        # Extract sequence using bedtools
        sequence=$(bedtools getfasta -fi "$fna_file" -bed <(echo -e "${scaffold}\t${start}\t${end}") 2>/dev/null | tail -n +2)

        # Reverse complement if it's the negative strand
        if [[ "$strand" == "-" ]]; then
            sequence=$(echo "$sequence" | rev | tr 'ACGTacgt' 'TGCAtgca')
        fi

        # Append the sequence to the specific gene type output file
        echo ">$geneid $gene_type $scaffold $start $end $strand $score" >> "$output_fasta"
        echo "$sequence" >> "$output_fasta"
    done

    echo "Sequences from $gene_file have been written to $output_fasta"
}

# Export the function to be available to srun
export -f process_gene_file
export genome_dir busco_dir output_dir

# Find only the files from 1st to 650th in sorted order
find "$busco_dir" -name "*_genes.txt" | sort | sed -n '1,650p' | xargs -P 64 -I {} srun --exclusive -n 1 bash -c 'process_gene_file "$1"' _ {}
