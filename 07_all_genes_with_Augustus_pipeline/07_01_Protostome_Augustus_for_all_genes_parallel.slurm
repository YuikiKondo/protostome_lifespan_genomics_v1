#!/bin/bash

## It is recommended to run around 100 species at a time as a single job. 
## Some species cannot be completed within 3 days if they have a large genome. In that case, split the genome file (species.fna) into several smaller files (e.g. 100 split), run AUGUSTUS for all the files, and combine the result files into a single output file (these processes are not described in this code).

#SBATCH --job-name=augustus_parallel
#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --mem=400G
#SBATCH --time=00-24:00:00
#SBATCH --output=07_01_Protostome_Augustus_for_all_genes_parallel.log
#SBATCH -p P112

# Load conda and activate AUGUSTUS environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate augustus_env

# Directories
busco_all="/home/y-kondo/protostome_lifespan_model_augustus/BUSCO_directory/BUSCO_all_species_copy"
augustus_species_dir="$CONDA_PREFIX/config/species"
genome_dir="/home/y-kondo/protostome_lifespan_model/ProtostomeGenomes_unzipped"
augustus_out_base="/home/y-kondo/protostome_lifespan_model_augustus/augustus_predictions/1_650"
completed_dir="/home/y-kondo/protostome_lifespan_model_augustus/augustus_predictions/completed_augustus"

mkdir -p "$augustus_out_base"

# === STEP 1: Get list of unprocessed species ===
all_targets=$(ls -d "$busco_all"/BUSCO_* | sed -n '1,650p')

species_to_process=()
for busco_dir in $all_targets; do
    genome_file="${busco_dir##*/BUSCO_}"
    species_model_name="${genome_file%.fna}"
    output_file="$completed_dir/${species_model_name}_augustus.gff"
    if [[ ! -f "$output_file" ]]; then
        species_to_process+=("$busco_dir")
    fi
done

# Output the species to process
echo "[INFO] Species to be processed by AUGUSTUS:"
printf "%s\n" "${species_to_process[@]}" | tee augustus_species_to_process.txt

# If none to process, exit early
if [[ ${#species_to_process[@]} -eq 0 ]]; then
    echo "[INFO] All species already processed. Exiting."
    exit 0
fi

# === STEP 2: Define function to run AUGUSTUS ===
run_augustus() {
    busco_dir="$1"
    species_tag=$(basename "$busco_dir")
    genome_file="${species_tag#BUSCO_}"
    genome_path="$genome_dir/$genome_file"
    species_model_name="${genome_file%.fna}"
    output_file="$completed_dir/${species_model_name}_augustus.gff"

    if [[ -f "$output_file" ]]; then
        echo "[SKIP] AUGUSTUS already completed for $species_model_name"
        return
    fi

    if [[ ! -f "$genome_path" ]]; then
        echo "[WARNING] Genome not found: $genome_path" >&2
        return
    fi

    model_dir="$busco_dir/run_metazoa_odb10/augustus_output/retraining_parameters"
    trained_model_path=$(find "$model_dir" -maxdepth 1 -type d -name "BUSCO_*")

    if [[ ! -d "$trained_model_path" ]]; then
        echo "[WARNING] No trained AUGUSTUS model found for $species_tag" >&2
        return
    fi

    target_model_dir="$augustus_species_dir/$species_model_name"
    mkdir -p "$target_model_dir"

    echo "[INFO] Copying AUGUSTUS model files to $target_model_dir"

    for f in "$trained_model_path"/*; do
        orig_filename=$(basename "$f")
        suffix=$(echo "$orig_filename" | sed "s/^.*\.fna_//")
        new_filename="${species_model_name}_${suffix}"
        echo "[DEBUG] Copying $orig_filename -> $new_filename"
        cp "$f" "$target_model_dir/$new_filename"
    done

    param_file="$target_model_dir/${species_model_name}_parameters.cfg"
    if [[ -f "$param_file" ]]; then
        echo "[INFO] Patching parameter file: $param_file"
        sed -i "s/BUSCO_.*\.fna_/${species_model_name}_/g" "$param_file"
    else
        echo "[WARNING] parameters.cfg file missing for $species_model_name" >&2
        return
    fi

    echo "[INFO] Running AUGUSTUS for $genome_file using model $species_model_name"
    augustus --species="$species_model_name" "$genome_path" > "$output_file"

    if [[ $? -ne 0 ]]; then
        echo "[ERROR] AUGUSTUS failed for $species_model_name" >&2
    else
        echo "[SUCCESS] Finished AUGUSTUS for $species_model_name"
    fi
}

export -f run_augustus
export augustus_species_dir genome_dir augustus_out_base completed_dir

# === STEP 3: Run AUGUSTUS in parallel only for species that need it ===
printf "%s\n" "${species_to_process[@]}" | xargs -P 56 -I {} srun --exclusive -n1 bash -c 'run_augustus "$@"' _ {}
