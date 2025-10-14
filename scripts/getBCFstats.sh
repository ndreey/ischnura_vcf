#!/bin/bash

#SBATCH --job-name vcf_stats
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

STATS_OUT=stats/vcf_stats.tsv


#########################################################################
#                  Summarize BCFTOOLS Stats
#########################################################################
# Get quick stats for the files
echo -e "$(date) [INFO]        Summarising stats from each VCF file: $STATS_OUT"
echo -e "VCF\tsamples\trecords\tinvariants\tSNPs\tMNPs\tindels\tothers\tmulti-sites\tmulti-SNPs" > $STATS_OUT

for FILE in $(ls vcfs/*.vcf.gz); do

    echo -e "$(date) [INFO]        Summarising $FILE"
    echo -e "$(basename $FILE)\t$(bcftools stats $FILE | grep -v "# SN," | grep -A 9 "# SN" | grep -v "#" | cut -f 4 | paste -sd$'\t')" >> $STATS_OUT

done

# Preview the file
echo -e ">> BCFTOOLS STATS DESCRIPTION << \n"
echo "number of records   .. number of data rows in the VCF"
echo "number of no-ALTs   .. reference-only sites, ALT is either '.' or identical to REF"
echo "number of SNPs      .. number of rows with a SNP"
echo "number of MNPs      .. number of rows with a MNP, such as CC>TT"
echo "number of indels    .. number of rows with an indel"
echo "number of others    .. number of rows with other type, for example a symbolic allele or a complex substitution, such as ACT>TCGA"
echo "number of multiallelic sites     .. number of rows with multiple alternate alleles"
echo "number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs"
echo "Note that rows containing multiple types will be counted multiple times, in each counter. For example, a row with a SNP and an indel increments both the SNP and the indel counter."
echo ""

# Print out the stats in a clean table
echo -e "\n >> BCFTOOLS STATS <<"
column -t "$STATS_OUT"



##########################################################################
##                                  VCF Stats
##########################################################################
#
#
#echo "$(date) [INFO]        Collecting site-level QC stats"
#
## RAW BIALLELIC
#if [ ! -f "$RAW_STATS_OUT" ]; then
#    echo "$(date) [INFO]        Writing raw biallelic site stats: $RAW_STATS_OUT"
#
#    echo -e "CHROM\tPOS\tQUAL\tDP\tQD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum\tAN\tAC\tAF" > "$RAW_STATS_OUT"
#
#    bcftools query -f '%CHROM\t%POS\t%QUAL\t%DP\t%QD\t%FS\t%SOR\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%AN\t%AC\t%AF\n' \
#        "$RAW_BIALLELIC" >> "$RAW_STATS_OUT"
#else
#    echo "$(date) [SKIP]        Raw biallelic site stats already exist: $RAW_STATS_OUT"
#fi
#
## Filtered biallelic
#if [ ! -f "$FILT_STATS_OUT" ]; then
#    echo "$(date) [INFO]        Writing raw biallelic site stats: $FILT_STATS_OUT"
#
#    echo -e "CHROM\tPOS\tQUAL\tDP\tQD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum\tAN\tAC\tAF" > "$FILT_STATS_OUT"
#
#    bcftools query -f '%CHROM\t%POS\t%QUAL\t%DP\t%QD\t%FS\t%SOR\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%AN\t%AC\t%AF\n' \
#        "$FILT_BIAL" >> "$FILT_STATS_OUT"
#else
#    echo "$(date) [SKIP]        Raw biallelic site stats already exist: $FILT_STATS_OUT"
#fi
#
## Get vcf stats
#for VCF_STAT in $RAW_BIALLELIC $FILT_BIAL; do
#
#    ID=$(basename $VCF_STAT .vcf.gz)
#    STAT_OUT=$STATSDIR/${PREFIX}-$ID
#    vcftools --gzvcf "$VCF_STAT" --site-mean-depth --out "$STAT_OUT"          # Get the per-site mean depth average across samples (.ldepth.mean)
#    vcftools --gzvcf "$VCF_STAT" --missing-site --out "$STAT_OUT"             # Get the proportion of missingness per site (.lmiss)
#    vcftools --gzvcf "$VCF_STAT" --depth --out "$STAT_OUT"                    # Get mean depth per sample (.idepth)
#    vcftools --gzvcf "$VCF_STAT" --missing-indv --out "$STAT_OUT"             # Get the proportion of missingness per sample (.imiss)
#    vcftools --gzvcf "$VCF_STAT" --SNPdensity 10000 --out "$STAT_OUT"         # Get SNP density across the genome (.snpden)
#done

echo -e "\n$(date) [FINISH]       Script Complete!\n"

