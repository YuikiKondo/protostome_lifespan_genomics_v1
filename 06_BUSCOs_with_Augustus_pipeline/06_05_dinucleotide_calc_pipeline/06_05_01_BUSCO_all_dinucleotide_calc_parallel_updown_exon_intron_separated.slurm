#!/bin/bash

#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --time=00-24:00:00
#SBATCH --output=/home/y-kondo/protostome_lifespan_model_augustus/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/dinucleotide_count/06_05_01_BUSCO_all_dinucleotide_calc_parallel_updown_exon_intron_separated.slurm.log
#SBATCH -p P128

echo "Job started at $(date)"
echo "Running on $(hostname)"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"

# Load Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate busco_env

# Run the Python script
python 06_05_01_BUSCO_all_dinucleotide_calc_parallel_updown_exon_intron_separated.py

echo "Job finished at $(date)"
