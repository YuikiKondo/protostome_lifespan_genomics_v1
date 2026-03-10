#!/bin/bash

#SBATCH --account=y-kondo
#SBATCH --mail-user=kondo-yuiki746@g.ecc.u-tokyo.ac.jp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1                       
#SBATCH --ntasks=1                     
#SBATCH --cpus-per-task=32               
#SBATCH --mem-per-cpu=100G
#SBATCH --time=00-24:00:00
#SBATCH --output=05_10_19_Rcode_PGLS_dinuc_raw_freq_10motifs_for_whole_genome_and_intergenic_with_genome_size_control_with_interaction_8Regions_using_nlme_gls_freq_as_response.slurm.log
#SBATCH -p P32

# Activate the Conda environment with R
source ~/miniconda3/bin/activate r_env

# Run the R script
Rscript 05_10_19_Rcode_PGLS_dinuc_raw_freq_10motifs_for_whole_genome_and_intergenic_with_genome_size_control_with_interaction_8Regions_using_nlme_gls_freq_as_response.r