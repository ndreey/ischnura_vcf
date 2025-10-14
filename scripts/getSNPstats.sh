#!/bin/bash

#SBATCH --job-name SNPstats
#SBATCH -A naiss2024-5-505
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=12GB
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

# Arguments
CPU=4
VCF=$1
PREFIX=$(basename $VCF .vcf.gz)
STATS_DIR=stats/$PREFIX
STATS_OUT=$STATS_DIR/raw_stats-$PREFIX.tsv

mkdir -p $STATS_DIR

echo "$(date) [INFO]        Writing raw biallelic site stats: $STATS_OUT"
# Write header line first
echo -e "CHROM\tPOS\tQUAL\tAN\tAC\tAF\tDP\tQD\tMQ\tFS\tSOR\tExcessHet\tInbreedingCoeff\tBaseQRankSum\tMQRankSum\tReadPosRankSum" > "$STATS_OUT"

# Query the info fields
bcftools query -f '%CHROM\t%POS\t%QUAL\t%INFO/AN\t%INFO/AC\t%INFO/AF\t%INFO/DP\t%INFO/QD\t%INFO/MQ\t%INFO/FS\t%INFO/SOR\t%INFO/ExcessHet\t%INFO/InbreedingCoeff\t%INFO/BaseQRankSum\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\n' \
    $VCF >> $STATS_OUT

# Get vcf stats
vcftools --gzvcf $VCF --missing-site --out $STATS_DIR/$PREFIX             # Get the proportion of missingness per site (.lmiss)
vcftools --gzvcf $VCF --depth --out $STATS_DIR/$PREFIX                    # Get mean depth per sample (.idepth)
vcftools --gzvcf $VCF --missing-indv --out $STATS_DIR/$PREFIX             # Get the proportion of missingness per sample (.imiss)
vcftools --gzvcf $VCF --freq2 --out $STATS_DIR/$PREFIX                     # Get the frequencies and MAF
echo -e "\n$(date) [FINISH]       Script Complete!\n"
