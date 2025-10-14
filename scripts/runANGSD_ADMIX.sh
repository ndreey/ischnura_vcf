#!/bin/bash

#SBATCH --job-name NGSadmix
#SBATCH -A naiss2024-5-505
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=56GB
#SBATCH -t 5:30:00
#SBATCH --output=slurm-logs/angsd/SLURM-%j.out
#SBATCH --error=slurm-logs/angsd/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date) [START]     Starting script execution"

# Load in modules
ml load angsd/0.940
ml load bcftools/1.20

# Arguments
VCF=$1
LD_KEEP=$2
SUFFIX=$3
ID=$(basename $VCF .vcf.gz)
BEAGLE=vcfs/$ID
CPU=8
ANGSD_DIR=04-ANGSD/$ID-$SUFFIX
LD_FILTER=$ANGSD_DIR/$ID.LD-prune-beagle.txt
LD_PRUNED_BEAGLE=$ANGSD_DIR/$ID.LDpruned.beagle.gz

mkdir -p $ANGSD_DIR

#echo "[INFO]    converting to beagle"
#angsd -vcf-pl $VCF -nThreads $CPU -doGlf 2 -doMajorMinor 1 -out $BEAGLE

echo "[INFO]    Convert PLINK LD prune file to beagle friendly format: $(basename $LD_KEEP)"
sed 's/:/_/g' $LD_KEEP > $LD_FILTER

echo "[INFO]    LD prune: $BEAGLE.beagle.gz"
{
  # keep header
  zcat "$BEAGLE.beagle.gz" | head -n1
  # grep exact marker IDs from the first column
  zcat "$BEAGLE.beagle.gz" | grep -F -f "$LD_FILTER"
} | bgzip > "$LD_PRUNED_BEAGLE"


echo "[INFO]    Running convergence"
NGSadmix \
    -likes $LD_PRUNED_BEAGLE \
    -K 2 \
    -o $ANGSD_DIR/$ID.K2 \
    -P $CPU

echo "$(date) [EXIT]     Complete"