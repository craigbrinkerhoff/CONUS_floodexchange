#!/bin/bash
#SBATCH --job-name='connectivity'
#SBATCH -c 10
#SBATCH -p day
#SBATCH --mail-type=END,FAIL
#SBATCi --mail-user=craig.brinkerhoff@yale.edu
#SBATCH --mem-per-cpu=10G #18G #30G
#SBATCH -t 24:00:00 
#SBATCH -o out_master.txt
#SBATCH -e err_master.txt

Rscript src/runTargets.R
