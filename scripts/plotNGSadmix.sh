#!/bin/bash

#SBATCH --job-name R_plot
#SBATCH -A naiss2024-5-505
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=12GB
#SBATCH -t 00:30:00
#SBATCH --output=slurm-logs/R/SLURM-%j.out
#SBATCH --error=slurm-logs/R/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date) [START]     Starting script execution"

# Load in the r-arena mamba environment
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate r-arena


FILE=nanala.txt
OUTDIR=test_plots

Rscript scripts/r-scripts/plotNGSadmix.R "$FILE" "$OUTDIR"
