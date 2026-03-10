#!/bin/bash

#SBATCH --account=username
#SBATCH --mail-user=your.email@domain.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=5  # Run 5 BUSCO jobs in parallel
#SBATCH --cpus-per-task=12  # Each job gets 10 CPUs
#SBATCH --mem=400G  # Adjust based on actual usage (was 800G)
#SBATCH --time=00-24:00:00
#SBATCH --output=/home/username/protostome_lifespan_model_augustus/BUSCO_directory/1_650/BUSCO_protostome_1_650_combined.log
#SBATCH -p P64

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# Source Conda initialization and activate augustus_env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate augustus_env

# Define paths
genome_dir="/home/username/protostome_lifespan_model/ProtostomeGenomes_unzipped"
busco_db="/home/username/protostome_lifespan_model_augustus/BUSCO_directory/busco_downloads/lineages/metazoa_odb10"
log_dir="/home/username/protostome_lifespan_model_augustus/BUSCO_directory/1_650/logs"
final_log="/home/username/protostome_lifespan_model_augustus/BUSCO_directory/1_650/BUSCO_protostome_1_650_combined.log"

mkdir -p "$log_dir"  # Ensure log directory exists

# Read genome filenames from the list into an array
genome_list="/home/username/protostome_lifespan_model/BUSCO_directory/genome_list_all_protostomes.txt"
mapfile -t genomes < <(sed -n '1,650p' "$genome_list")	# Process files from 1st to 650th.(If you want to process fewer files, change '1,650p'.)

# Function to run BUSCO
run_busco() {
    genome_file="${genome_dir}/$1"
    log_file="${log_dir}/$(basename $1).log"
    
    echo "=============================================" > "$log_file"
    echo "Processing file: ${genome_file}" >> "$log_file"
    echo "Running BUSCO on: ${genome_file}" >> "$log_file"

    # Run BUSCO and write to a separate log file
    busco -i "${genome_file}" -l "${busco_db}" -m genome --long -c 12 --augustus >> "$log_file" 2>&1

    # Append to the final combined log after completion
    cat "$log_file" >> "$final_log"
    echo "" >> "$final_log"  # Add an empty line for clarity

    # Optional: Remove individual log file after merging (uncomment below if needed)
    # rm "$log_file"
}

export -f run_busco
export genome_dir busco_db log_dir final_log

# Run multiple BUSCO jobs in parallel using Slurm's built-in job array mechanism
for genome in "${genomes[@]}"; do
    srun --exclusive -n1 --cpus-per-task=12 bash -c "run_busco '$genome'" &
done

wait  # Wait for all jobs to complete

echo "All BUSCO analyses completed. Final log saved at: $final_log"
