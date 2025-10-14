#!/bin/bash

#SBATCH --job-name convergence
#SBATCH -A naiss2024-5-505
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=46GB
#SBATCH -t 24:30:00
#SBATCH --output=slurm-logs/convergence/SLURM-%j.out
#SBATCH --error=slurm-logs/convergence/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date) [START]     Starting script execution"

# Load in modules
ml load bioinfo-tools
ml load plink/2.00a5.14
ml load ADMIXTURE/1.3.0

# Arguments
CPU=12
VCF=$1
ID=$(basename $VCF .vcf.gz)
WINDOW=$2
STEP=$3
R=$4
PREFIX=$ID-$WINDOW-$STEP-$R
LD_OUT=001-pruneLD/$PREFIX
PCA_OUT=002-PCA/$PREFIX
ADMIX_BURN=0004-ADMIX_CONVERGENCE/$PREFIX
WORK=$(pwd)

mkdir -p $LD_OUT $PCA_OUT $ADMIX_BURN


echo "$(date) [INFO]     LD pruning"
plink2 \
    --vcf $VCF \
    --threads $CPU \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids "@:#" \
    --indep-pairwise $WINDOW $STEP $R \
    --out $LD_OUT/$PREFIX

echo "$(date) [EXEC]     Run PCA"
plink2 \
    --vcf $VCF \
    --threads $CPU \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids "@:#" \
    --extract $LD_OUT/$PREFIX.prune.in \
    --make-bed \
    --pca \
    --out $PCA_OUT/$PREFIX


echo "$(date) [Info]     Setting up for ADMIXTURE"
# ADMIXTURE does not accept chromosome names that are not human chromosomes.
# Change name of original
mv $PCA_OUT/$PREFIX.bim $PCA_OUT/$PREFIX.human-chr.bim

# Create a new one with first column set to 0
cat $PCA_OUT/$PREFIX.human-chr.bim | \
    awk '{$1="0"; OFS="\t"; print $0}' > $PCA_OUT/$PREFIX.bim


echo "$(date) [Info]     Init convergence burn (max 100 runs)"
# EXPLANATION
#
#
#
#

# Folders to put individual runs.
mkdir -p "$ADMIX_BURN"/{logs,Qs,Ps}
cd "$ADMIX_BURN"

for K in {1..3}; do
    # Arguments
    MAX_RUNS=100
    THRESH=2.0
    runs=0
    BURN_TSV=$K-burn-summary.tsv

    # Start an empty .tsv
    : > $BURN_TSV
    
    while [ $runs -lt $MAX_RUNS ]; do
        runs=$((runs+1))
        log=$PREFIX-log-K_$K-run_$runs.out

        echo "[K=$K run=$runs] starting..."
        admixture --cv -j8 -s time ../../$PCA_OUT/$PREFIX.bed $K > $log

        # keep Q/P for this run so they don't overwrite each other
        mv "${PREFIX}.${K}.Q" "Qs/${PREFIX}.K${K}.r${runs}.Q"
        mv "${PREFIX}.${K}.P" "Ps/${PREFIX}.K${K}.r${runs}.P"

        # Grab the loglikelihood value
        logLH=$(cat $log | grep -v "Elapsed" | grep -Eo 'Loglikelihood:[[:space:]]+-?[0-9.]+([eE][-+]?[0-9]+)?' | cut -f 2 -d " ")
        
        # Write to tsv.
        echo -e "$K\t$runs\t$logLH\t$log" >> $BURN_TSV

        if [ $runs -ge 3 ]; then
            # extract top-3 LLs for this K (highest = least negative)
            top3=$(awk -v k=$K '$1==k {print $3}' $BURN_TSV | sort -g | tail -n 3)
            best=$(echo "$top3" | tail -n 1)
            third=$(echo "$top3" | head -n 1)
            diff=$(echo "$best - $third" | bc -l)

            echo "[K=$K run=$runs] best=$best third=$third diff=$diff"

            check=$(echo "$diff <= $THRESH" | bc -l)
            if [ "$check" -eq 1 ]; then
                echo "[K=$K] converged after $runs runs (diff ≤ $THRESH)."
                break
            fi
        fi
    done

    if [ $runs -eq $MAX_RUNS ]; then
        echo "[K=$K] reached MAX_RUNS=$MAX_RUNS without convergence."
    fi
done

echo "$(date) [COMPLETE]  Done!"
