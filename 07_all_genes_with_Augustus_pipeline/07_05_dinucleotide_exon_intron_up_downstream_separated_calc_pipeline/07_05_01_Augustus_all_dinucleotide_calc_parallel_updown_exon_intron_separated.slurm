#!/bin/bash

#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=400G
#SBATCH --time=00-24:00:00
#SBATCH --output=07_05_01_Augustus_all_dinucleotide_calc_parallel_updown_exon_intron_separated.slurm.log
#SBATCH -p P112

echo "Job started at $(date)"
echo "Running on $(hostname)"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"

# Load Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate busco_env

# Run the Python script
python 07_05_01_Augustus_all_dinucleotide_calc_parallel_updown_exon_intron_separated.py

echo "Job finished at $(date)"
