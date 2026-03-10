#!/bin/bash

#SBATCH --account=your_account
#SBATCH --mail-user= your.email@domain.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=20G
#SBATCH --time=00-08:00:00
#SBATCH --output=%x.%j.log
#SBATCH -p P64

## create variables for directories and files
WD=/home/username/protostome_lifespan_model                ## working directory
cd $WD

bash decompress_fna.sh

## END SCRIPT
