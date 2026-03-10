#!/bin/bash

#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1                       
#SBATCH --ntasks=1                     
#SBATCH --cpus-per-task=40               
#SBATCH --mem-per-cpu=20G
#SBATCH --time=00-24:00:00
#SBATCH --output=/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/TSS_Github_check/09_03_sliding_window_all_dinucleotide_counts_rates_parallel_precise_calc.R.log
#SBATCH -p P128

# Activate the Conda environment with R
source ~/miniconda3/bin/activate r_env

# Run the R script
Rscript /home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/TSS_Github_check/09_03_sliding_window_all_dinucleotide_counts_rates_parallel_precise_calc.R