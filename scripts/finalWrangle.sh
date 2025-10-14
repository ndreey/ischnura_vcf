#!/bin/bash

#SBATCH --job-name finalWrangle
#SBATCH -A naiss2024-5-505
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=56GB
#SBATCH -t 16:30:00
#SBATCH --output=slurm-logs/plink/SLURM-%j.out
#SBATCH --error=slurm-logs/plink/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date) [START]     Starting script execution"

# Load in modules
ml load bcftools/1.20

# Arguments
CPU=8
VCF=vcfs/andre.noX13.filt.bial.vcf.gz
VCF_OUT=vcfs/andre.noX13.final.bial.vcf.gz
PREFIX=$(basename $VCF_OUT .vcf.gz)

echo "$(date) [INFO]     Filtering on site missingngess and MAF"
# Exclude records with 20% or more missigness and MAF below 0.02.
bcftools filter --threads $CPU -e 'F_MISSING > 0.2 || MAF <= 0.03' $VCF \
| bcftools view --threads $CPU -f PASS -O z -o $VCF_OUT 

echo "$(date) [INFO]     Index vcf"
tabix $VCF_OUT

echo "$(date) [EXEC]     Get overall stats"
sbatch scripts/getBCFstats.sh
sbatch scripts/getSNPstats.sh $VCF_OUT

echo "$(date) [COMPLETE]     Done!"