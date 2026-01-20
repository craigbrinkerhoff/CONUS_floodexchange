#!/bin/bash
#SBATCH --job-name='connectivity'
#SBATCH -c 10 #10 # Number of Cores per Task
#SBATCH -p week
#SBATCH --mail-type=END,FAIL
#SBATCi --mail-user=craig.brinkerhoff@yale.edu
#SBATCH --mem-per-cpu=10G #18G #30G #Requested memory
#SBATCH -t 72:00:00  # Job time limit
#SBATCH -o out_master.txt  # %j = job ID
#SBATCH -e err_master.txt

Rscript src/runTargets.R
