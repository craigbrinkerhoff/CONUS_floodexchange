#!/bin/bash
#SBATCH --job-name='connectivity'
#SBATCH -c 20  # Number of Cores per Task
#SBATCH -p day
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=craig.brinkerhoff@yale.edu
#SBATCH --mem=50000 #Requested memory
#SBATCH -t 24:00:00  # Job time limit
#SBATCH -o out_master.txt  # %j = job ID
#SBATCH -e err_master.txt

Rscript src/runTargets.R
