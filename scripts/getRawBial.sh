#!/bin/bash

#SBATCH --job-name bcftools
#SBATCH -A naiss2024-5-505
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=36GB
#SBATCH -t 04:30:00
#SBATCH --output=slurm-logs/wrangle/SLURM-%j.out
#SBATCH --error=slurm-logs/wrangle/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date) [START]     Starting script execution"

# Laad modules
ml load vcftools/0.1.16
ml load bcftools/1.20

###### Define variables
# Raw Variants to extract stats from
VCF_RAW=../11_vcf_wrangling/vcf-final/SwD_VCF.variants.vcf.gz
# Output
RAW_BIAL=vcfs/SwD_VCF.raw.bial.vcf.gz
CPU=8

# Create comma separated string of regions to exclude.
EXCLUDE=$(paste -sd "," doc/excluded_chr.txt)

# Get biallelic snps of all regions besides excluded regions.
# -m2: minimum number of alleles = 2
# -M2: maximum number of alleles = 2
# -v snps: only select snps.
echo "$(date) [INFO]     Subset vcf"
bcftools view \
    --threads $CPU \
    -r $EXCLUDE \
    -m2 \
    -M2 \
    -v snps \
    -O z \
    -o $RAW_BIAL \
    $VCF_RAW

echo "$(date) [INFO]     Indexing vcf"
tabix --threads $CPU $RAW_BIAL

echo "$(date) [FINISH]     Done"
