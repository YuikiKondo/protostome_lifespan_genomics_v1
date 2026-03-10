#!/bin/bash

#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=400G
#SBATCH --time=00-24:00:00
#SBATCH --output=/home/y-kondo/protostome_lifespan_model/Genome_dinucleotide_count/04_01_01_genome_all_dinucleotide_calc_parallel4.slurm.log
#SBATCH -p P112

echo "Job started at $(date)"
echo "Running on $(hostname)"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"

# Load Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate busco_env

# Run the Python script
python 04_01_01_genome_all_dinucleotide_calc_parallel4.py

echo "Job finished at $(date)"