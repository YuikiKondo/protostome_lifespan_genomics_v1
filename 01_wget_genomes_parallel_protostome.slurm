#!/bin/bash

#SBATCH --account=your_account
#SBATCH --mail-user=your.email@domain.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00-08:00:00
#SBATCH --output=%x.%j.log
#SBATCH -p SL_A

## SET UP ----

## create variables for directories and files
export WD=/home/username/protostome_lifespan_model	## working directory
export DD=$WD/dataFiles			                ## data directory
export ODT=$WD/scratch/tempProtostomeGenomes	## temporary output directory
export OD=$PWD/ProtostomeGenomes			## output directory

## make the output directories
mkdir -p $ODT
mkdir -p $OD

## move to scratch
cd scratch

## Add directories to PATH
export PATH=$HOME/.local/bin:$HOME/go/bin:$PATH

echo "Current PATH: $PATH"
which parallel
which pget

## set the path to the TSV file
TSV_FILE="$DD/00.00_genbank_best_protostome_genomes_dl.csv"


## run 60 jobs at a time in parallel
cat "$TSV_FILE" | parallel --header : --jobs 60 --colsep ',' '

  ## create variables from common names
  organism_name={1}
  assembly_id={2}
  genome_ftp_path={3}
  checksums_ftp_path={4}

  genome_file_name=$(echo $genome_ftp_path | cut -d "/" -f 11)
  echo ${genome_file_name}

  checksums_file_name=$(echo $checksums_ftp_path | \
  sed "s#/md5checksums.txt#_md5checksums.txt#g" |
  cut -d "/" -f 10)

  if [ -e $OD/$genome_file_name ]; then
	echo $organism_name" genome ("$genome_file_name") is already downloaded"
  else
        echo $organism_name" genome (assembly "$assembly_id") and checksums ** downloading now **"
        pget -p 4 $genome_ftp_path -o $ODT/$genome_file_name
        ## grep the checksums value for the genome only
        wget $checksums_ftp_path -O - | grep $genome_file_name > $ODT/$checksums_file_name
  fi
'

## rsync folder to datasets
rsync -avh "$ODT/" "$OD" && \
echo "Rsync completed successfully" || \
echo "Rsync failed, check the log file for more information"

## remove folder from scratch
rm -r "$ODT" && echo "Source directory deleted"
echo "Job complete"

## END SCRIPT