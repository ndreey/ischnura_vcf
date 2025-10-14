#!/bin/bash

#SBATCH --job-name bcftools
#SBATCH -A naiss2024-5-505
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=36GB
#SBATCH -t 06:30:00
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
VCF_RAW=vcfs/SwD_VCF.raw.bial.vcf.gz

# Filtered but not biallelic
VCF_IN=../11_vcf_wrangling/vcf-final/SwD_VCF.filt.no-miss.vcf.gz

# Temporary that holds all the snps we want
VCF_TMP=vcfs/SwD_VCF.filt.no-miss.no-X13.bial.tmp.vcf.gz

# Biallelic without CHR X or 13, with base filtering.
VCF_OUT=vcfs/SwD_VCF.filt.no-miss.no-X13.bial.vcf.gz

# Positions to restrict to
KEEP_POS=doc/standard_filter_pos.txt

CPU=8

# Create comma separated string of regions to exclude.
EXCLUDE=$(paste -sd "," doc/excluded_chr.txt)

# Get biallelic snps of all regions besides excluded regions.
# -m2: minimum number of alleles = 2
# -M2: maximum number of alleles = 2
# -v snps: only select snps.
#echo "$(date) [INFO]     Subset vcf"
#bcftools view \
#   --threads $CPU \
#    -t $EXCLUDE \
#    -m2 \
#    -M2 \
#    -v snps \
#    -O z \
#    -o $VCF_TMP \
#    $VCF_IN

echo "$(date) [INFO]     Indexing vcf"
#tabix --threads $CPU $VCF_TMP
bcftools index --threads $CPU $VCF_TMP

echo "$(date) [INFO]     Extract chr and pos"
bcftools query -f '%CHROM\t%POS\n' $VCF_TMP > $KEEP_POS

echo "$(date) [INFO]     Indexing raw vcf"
bcftools index --threads $CPU $VCF_RAW

echo "$(date) [INFO]     Extract bial snps from raw vcf that holds stats"
# Extract from raw variants to retain the stats info.
bcftools view --threads $CPU -T $KEEP_POS -O z -o $VCF_OUT $VCF_RAW

echo "$(date) [INFO]     Indexing final vcf: $VCF_OUT"
tabix --threads $CPU $VCF_OUT
bcftools index --threads $CPU $VCF_OUT

echo "$(date) [FINISH]     Done!"