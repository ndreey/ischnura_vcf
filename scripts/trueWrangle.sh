#!/bin/bash

#SBATCH --job-name bcftools
#SBATCH -A naiss2024-5-505
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=56GB
#SBATCH -t 16:30:00
#SBATCH --output=slurm-logs/filter/SLURM-%j.out
#SBATCH --error=slurm-logs/filter/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date) [START]     Starting script execution"

ml load bcftools/1.20

# Arguments
CPU=8
VARIANTS_VCF=../11_vcf_wrangling/vcf-final/SwD_VCF.variants.vcf.gz
VCF_BIAL_HARD=vcfs/andre.noX13.hard.bial.vcf.gz
VCF_FILTERED=vcfs/andre.noX13.filt.bial.vcf.gz

echo "$(date) [INFO]     Create index"
# .tbi index exists, but lets also create .csi (might be quicker lookup??)

# Create comma separated string of regions to exclude.
EXCLUDE=$(paste -sd "," doc/excluded_chr.txt)

echo "$(date) [INFO]     Hard filters + SnpGap + biallelic SNPs, excluding listed contigs"
# Apply GATK hard filtering, remove monomorphic ALT sites, and remove snps in 10bp proximity of indels
# Only keep biallelic snps that passed the filter in the regions besides chr X and 13
bcftools filter --threads $CPU --SnpGap 10 -e 'AC==0 || AC==AN || FS>60.0 || SOR>3 || MQ<40 || MQRankSum<-10.0 || MQRankSum>6.0 || QD<2.0 || ReadPosRankSum<-8.0' $VARIANTS_VCF \
| bcftools view --threads $CPU -t ^$EXCLUDE -m2 -M2 -v snps -f PASS -O z -o $VCF_BIAL_HARD

echo "$(date) [INFO]     Index the VCF!"
# Index the file
tabix $VCF_BIAL_HARD

echo "$(date) [INFO]     Apply individual genotype filtering to missing if below thresholds"
# Lets do further quality filtering
# Replace individual GT with missing if depth <3 or quality is < 20
bcftools filter --threads $CPU -S . -e 'FMT/DP<3 | FMT/GQ<20' -O z -o $VCF_FILTERED $VCF_BIAL_HARD

echo "$(date) [INFO]        Index final VCF!"
tabix $VCF_FILTERED

echo "$(date) [FINISH]     Filtering done!"