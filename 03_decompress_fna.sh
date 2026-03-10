#!/bin/bash

# Directory with compressed .fna.gz files
input_dir="/home/username/protostome_lifespan_model/ProtostomeGenomes/"
# Directory to save the decompressed .fna files
output_dir="/home/username/protostome_lifespan_model/ProtostomeGenomes_unzipped"

# Create the output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop through all .fna.gz files in the input directory
for gz_file in "${input_dir}"/*.fna.gz; do
    # Get the base name of the file (removing the .gz extension)
    base_name=$(basename "${gz_file}" .gz)
    
    # Check if the decompressed file already exists
    if [ -f "${output_dir}/${base_name}" ]; then
        echo "Skipping ${base_name}, already exists in ${output_dir}"
    else
        echo "Decompressing ${gz_file} to ${output_dir}/${base_name}"
        gunzip -c "${gz_file}" > "${output_dir}/${base_name}"
    fi
done

echo "Decompression completed. All missing files have been processed."
