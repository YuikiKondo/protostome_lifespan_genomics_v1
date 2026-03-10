#!/bin/bash

#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1                       
#SBATCH --ntasks=36                     
#SBATCH --cpus-per-task=1               
#SBATCH --mem-per-cpu=10G
#SBATCH --time=00-24:00:00
#SBATCH --output=/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/TSS_Github_check/09_01_Extract_busco_TSSsequences_updown_exon_intron_separted_parallel.slurm.log
#SBATCH -p P128

export PATH=$HOME/miniconda3/bin:$PATH
source activate busco_env


# representative_species directory includes three BUSCO result files (i.e. complete, duplicated, fragmented) of 12 species (shortest and longest-lived species in 6 phyla)
busco_dir="/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/TSS_Github_check/representative_species/"
genome_dir="/home/y-kondo/protostome_lifespan_model/ProtostomeGenomes_unzipped/"
output_dir="/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/TSS_Github_check/representative_species_TSS_coordinates/"

mkdir -p "$output_dir"

process_gene_file() {
    gene_file="$1"
    gene_filename=$(basename "$gene_file")

    echo "Processing gene file: $gene_file"

    species=$(basename "$gene_file" | sed -E 's/^BUSCO_//' | sed -E 's/_(complete|fragmented|duplicated)_genes\.txt$//')
    if [[ "$gene_file" =~ (complete|fragmented|duplicated) ]]; then
        gene_type="${BASH_REMATCH[1]}"
    else
        echo "Error: Gene type not found in file name"
        return
    fi

    output_fasta="${output_dir}/${species}_${gene_type}_TSS_regions.fasta"
    fna_file="${genome_dir}/${species}"

    if [[ ! -f "$fna_file" ]]; then
        echo "Error: Genomic file for species $species does not exist!"
        return
    fi

    > "$output_fasta"

    tail -n +2 "$gene_file" | while read -r line; do
        region_type=$(echo "$line" | awk '{print $1}')
        scaffold=$(echo "$line" | awk '{print $2}')
        start=$(echo "$line" | awk '{print $3}')
        end=$(echo "$line" | awk '{print $4}')
        strand=$(echo "$line" | awk '{print $5}')
        score=$(echo "$line" | awk '{print $6}')

        # Only process "_1_exon" lines to get TSS regions
        if [[ "$region_type" =~ "_1_exon" ]]; then
            if [[ "$strand" == "+" ]]; then
                tss_start=$((start - 3000))
                tss_end=$((start + 3000))
            else
                tss_start=$((end - 3000))
                tss_end=$((end + 3000))
            fi

            # Extract the sequence for the TSS region
            tss_sequence=$(bedtools getfasta -fi "$fna_file" -bed <(echo -e "${scaffold}\t${tss_start}\t${tss_end}") 2>/dev/null | tail -n +2)

            # Reverse complement if it's the negative strand
            if [[ "$strand" == "-" ]]; then
                tss_sequence=$(echo "$tss_sequence" | rev | tr 'ACGTacgt' 'TGCAtgca')
            fi

            echo ">$region_type TSS_region ${scaffold}:${tss_start}-${tss_end} $strand $score" >> "$output_fasta"
            echo "$tss_sequence" >> "$output_fasta"
        fi
    done

    echo "TSS sequences from $gene_file have been written to $output_fasta"
}

export -f process_gene_file
export genome_dir busco_dir output_dir

find "$busco_dir" -name "*_genes.txt" | sort | sed -n '1,36p' | parallel -j 36 process_gene_file
