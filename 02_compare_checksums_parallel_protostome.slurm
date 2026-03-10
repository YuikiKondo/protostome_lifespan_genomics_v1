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

## create variables for directories and files
WD=/home/username/protostome_lifespan_model                ## working directory
DD=$WD/ProtostomeGenomes              ## data directory

## load module
export PATH=$HOME/.local/bin:$PATH

## move to dir
cd $DD

## compare checksums
ls *.fna.gz | parallel --jobs 60 '
  GENOME={}
  ID=$(basename $GENOME _genomic.fna.gz)
  TXT=$(echo $ID"_md5checksums.txt")
  echo "ID: " $ID
  CS1=$(md5sum $GENOME | awk "{ print \$1 }")
  echo $CS1
  CS2=$(awk "{ print \$1 }" < $TXT)
  echo $CS2
  if [ "$CS1" == "$CS2" ]
  then
    echo "all good"
  else
    echo "check"
  fi
' | tee $WD/02_compare_checksums_results.txt

grep -B 3 'check' $WD/02_compare_checksums_results.txt > "$WD/02_to_check.$(date +%Y%m%d).txt"

rm $WD/02_compare_checksums_results.txt

## END SCRIPT
