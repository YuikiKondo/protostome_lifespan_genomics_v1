#!/bin/bash

# Set paths
input_dir="/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/busco_coordinates_including_dup"
csv_file="/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/TSS_Github_check/Assembly_species_match_revised.csv"
output_dir="/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/TSS_Github_check/representative_species"

# Make sure output directory exists
mkdir -p "$output_dir"

# Skip the header and loop through assembly IDs
tail -n +2 "$csv_file" | cut -d',' -f1 | while read -r assembly; do
    # Use wildcard matching to find all three types of BUSCO files
    for status in complete duplicated fragmented; do
        src_file="$input_dir/BUSCO_${assembly}_genomic.fna_${status}_genes.txt"
        if [[ -f "$src_file" ]]; then
            cp "$src_file" "$output_dir/"
            echo "Moved: $src_file"
        else
            echo "Not found: $src_file"
        fi
    done
done
