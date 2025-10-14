#!/bin/bash

#SBATCH --job-name R_plot
#SBATCH -A naiss2024-5-505
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=12GB
#SBATCH -t 02:30:00
#SBATCH --output=slurm-logs/R/SLURM-%j.out
#SBATCH --error=slurm-logs/R/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date) [START]     Starting script execution"

# Load in the r-arena mamba environment
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate r-arena

PREFIX=SwD_VCF.noX13.final.bial
DIR=stats/$PREFIX
OUT_DIR=$DIR/plots
EIGEN_VAL=002-PCA/SwD_VCF.noX13.final.bial.eigenval
EIGEN_VEC=002-PCA/SwD_VCF.noX13.final.bial.eigenvec

mkdir -p $DIR $OUT_DIR

echo "$(date) [Info]     Plotting PCA"
Rscript scripts/r-scripts/plotPCA.R \
    $EIGEN_VAL \
    $EIGEN_VEC \
    $OUT_DIR


echo "$(date) [DONE]       Complete!"