#!/bin/bash
#SBATCH --job-name='connectivity'
#SBATCH -c 5 #9
#SBATCH -p day
#SBATCH --mail-type=END,FAIL
#SBATCi --mail-user=craig.brinkerhoff@yale.edu
#SBATCH --mem=48G #18G #30G
#SBATCH -t 24:00:00 
#SBATCH -o out_master.txt
#SBATCH -e err_master.txt

Rscript src/runTargets.R
