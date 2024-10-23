#!/bin/bash
#SBATCH --job-name='test_conn'
#SBATCH -c 20  # Number of Cores per Task
#SBATCH -p ceewater_cjgleason-cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cbrinkerhoff@umass.edu
#SBATCH --mem=50000 #Requested memory
#SBATCH -t 9:00:00  # Job time limit
#SBATCH -o out_master.txt  # %j = job ID
#SBATCH -e err_master.txt

Rscript src/runTargets.R
