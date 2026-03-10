#!/bin/bash

#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=00-24:00:00
#SBATCH --output=07_03_02_augustus_extract_gene_coordinates_-2kb_to_+1kb_front_all_duplicated.slurm.log
#SBATCH -p P64

echo "Job started at $(date)"
echo "Running on $(hostname)"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"

# Load Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate augustus_env

# Run the Python script
python 07_03_02_augustus_extract_gene_coordinates_-2kb_to_+1kb_front_all_duplicated.py

echo "Job finished at $(date)"
