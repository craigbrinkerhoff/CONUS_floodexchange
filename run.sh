#!/bin/bash
#SBATCH --job-name='connectivity'
#SBATCH -c 8 # Number of Cores per Task
#SBATCH -p day
#SBATCH --mail-type=END,FAIL
#SBATCi --mail-user=craig.brinkerhoff@yale.edu
#SBATCH --mem=90G #Requested memory
#SBATCH -t 24:00:00  # Job time limit
#SBATCH -o out_master.txt  # %j = job ID
#SBATCH -e err_master.txt

Rscript src/runTargets.R
