#!/bin/bash
#SBATCH --job-name=NGSadmix
#SBATCH -A naiss2024-5-505
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=56GB
#SBATCH -t 15:00:00
#SBATCH --output=slurm-logs/angsd/SLURM-%j.out
#SBATCH --error=slurm-logs/angsd/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

#set -euo pipefail

# ----------------------------- #
# Inputs
# ----------------------------- #
VCF=$1      # gzipped VCF
WINDOW=$2   # e.g. 500
STEP=$3     # e.g. 50
R=$4        # e.g. 0.05

CPU=8
KLIST="1 2 3"           # adjust if you want more K's
MAX_RUNS=100
THRESH=2.0              # LL units: best - third <= THRESH => converged

ID=$(basename "$VCF" .vcf.gz)
TAG="${ID}-${WINDOW}-${STEP}-${R}"

# Layout
PRUNE_DIR="01-pruneLD/${TAG}"
BEAGLE_DIR="beagles"
WORK_DIR="04-ANGSD/${TAG}"
LOG_DIR="${WORK_DIR}/logs"
RUN_DIR="${WORK_DIR}/runs"
BEST_DIR="${WORK_DIR}/best"
SUM_DIR="${WORK_DIR}/summaries"

mkdir -p "$PRUNE_DIR" "$BEAGLE_DIR" "$WORK_DIR" "$LOG_DIR" "$RUN_DIR" "$BEST_DIR" "$SUM_DIR" slurm-logs/angsd

PLINK_PREFIX="${PRUNE_DIR}/${TAG}"
BEAGLE="${BEAGLE_DIR}/${ID}.beagle.gz"
LD_FILTER="${PRUNE_DIR}/${TAG}.LD-prune-beagle.txt"
LD_PRUNED="${BEAGLE_DIR}/${TAG}.LDpruned.beagle.gz"
SAMPLES_TXT="${WORK_DIR}/${ID}.samples.txt"
MASTER_LOG="${SUM_DIR}/convergence.log"

echo "$(date) [START] Pipeline for $ID ($WINDOW/$STEP/$R)"

# ----------------------------- #
# Modules
# ----------------------------- #
ml load plink/2.00a5.14
ml load angsd/0.940
ml load bcftools/1.20

# ----------------------------- #
# Step 0: Samples file (once)
# ----------------------------- #
if [[ ! -s "$SAMPLES_TXT" || "$SAMPLES_TXT" -ot "$VCF" ]]; then
  echo "[INFO] Extracting sample IDs from VCF -> $SAMPLES_TXT"
  bcftools query -l "$VCF" > "$SAMPLES_TXT"
fi

# ----------------------------- #
# Step 1: PLINK LD prune (idempotent)
# ----------------------------- #
if [[ -s "${PLINK_PREFIX}.prune.in" && -s "${PLINK_PREFIX}.prune.out" ]]; then
  echo "[INFO] Reusing existing PLINK prune files in $PRUNE_DIR"
else
  echo "[INFO] Running PLINK LD prune -> $PRUNE_DIR"
  plink2 \
    --vcf "$VCF" \
    --threads "$CPU" \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids "@:#" \
    --indep-pairwise "$WINDOW" "$STEP" "$R" \
    --out "$PLINK_PREFIX"
fi

# Make Beagle-style keep list (underscore IDs)
if [[ ! -s "$LD_FILTER" || "$LD_FILTER" -ot "${PLINK_PREFIX}.prune.in" ]]; then
  echo "[INFO] Converting prune list to Beagle IDs -> $LD_FILTER"
  sed 's/:/_/g' "${PLINK_PREFIX}.prune.in" > "$LD_FILTER"
fi

# ----------------------------- #
# Step 2: Beagle (idempotent)
# ----------------------------- #
if [[ -s "$BEAGLE" ]]; then
  echo "[INFO] Reusing existing Beagle: $BEAGLE"
else
  echo "[INFO] Creating Beagle with ANGSD -> $BEAGLE"
  angsd -vcf-pl "$VCF" -nThreads "$CPU" -doGlf 2 -doMajorMinor 1 -out "${BEAGLE_DIR}/${ID}"
fi

# ----------------------------- #
# Step 3: LD-pruned Beagle (existence-only)
# ----------------------------- #
if [[ -f "$LD_PRUNED" ]]; then
  echo "[INFO] Reusing existing LD-pruned Beagle: $LD_PRUNED"
else
  echo "[INFO] LD prune Beagle -> $LD_PRUNED"
  {
    zcat "$BEAGLE" | head -n1
    zcat "$BEAGLE" | grep -F -f "$LD_FILTER"
  } | bgzip > "$LD_PRUNED"
fi


# ----------------------------- #
# Step 4: NGSadmix convergence loop (pure bash)
# ----------------------------- #

# master log header (append-safe)
if [[ ! -s "$MASTER_LOG" ]] || ! grep -q '^datetime' "$MASTER_LOG"; then
  echo -e "datetime\tID\tK\trun\tbestLL\tthirdLL\tdiffLL\tstatus" >> "$MASTER_LOG"
fi

echo "does it work?"

for K in $KLIST; do
  echo "[INFO] ==== NGSadmix K=$K ===="
  BURN_TSV="${SUM_DIR}/K${K}-burn-summary.tsv"
  echo -e "K\trun\tLL\tlog\tqopt\tfopt" > "$BURN_TSV"

  runs=0
  converged=0

  while [[ $runs -lt $MAX_RUNS ]]; do
    runs=$((runs+1))
    TAG_RUN="${ID}.K${K}.r${runs}"
    OUTPFX="${RUN_DIR}/${TAG_RUN}"
    RLOG="${LOG_DIR}/${TAG_RUN}.log"

    echo "[K=$K run=$runs] start"
    NGSadmix \
      -likes "$LD_PRUNED" \
      -K "$K" \
      -outfiles "$OUTPFX" \
      -P "$CPU" \
      -minMaf 0 \
      -misTol 0.05 \
      > "$RLOG" 2>&1

    # Parse final log-likelihood
    LL=$(grep -o 'best like=[^ ]*' "$RLOG" | tail -n1 | cut -d= -f2)
    if [[ -z "${LL:-}" ]]; then
      # fallback to last iter[...] like is=...
      LL=$(grep -o 'iter\[[0-9]\+\] like is=[^ ]*' "$RLOG" | tail -n1 | awk '{print $4}')
    fi

    QFILE="${OUTPFX}.qopt"
    FFILE="${OUTPFX}.fopt.gz"
    echo -e "${K}\t${runs}\t${LL}\t${RLOG}\t${QFILE}\t${FFILE}" >> "$BURN_TSV"

    # Convergence check once we have ≥3 runs: top1 vs top3 ≤ THRESH
    nrows=$(awk 'END{print NR-1}' "$BURN_TSV")
    if [[ "$nrows" -ge 3 ]]; then
      # best = largest LL; third = smallest among the top3 LLs
      best=$(awk 'NR>1{print $3}' "$BURN_TSV" | sort -g | tail -n 1)
      third=$(awk 'NR>1{print $3}' "$BURN_TSV" | sort -g | tail -n 3 | head -n 1)
      diff=$(echo "$best - $third" | bc -l)

      if echo "$diff <= $THRESH" | bc -l | grep -q '^1$'; then
        status="converged"
        converged=1
      else
        status="continue"
      fi
      echo -e "$(date +'%F %T')\t$ID\t$K\t$runs\t$best\t$third\t$diff\t$status" | tee -a "$MASTER_LOG"

      if [[ "$converged" -eq 1 ]]; then
        echo "[K=$K] Converged after $runs runs (diff ≤ $THRESH)."
        break
      fi
    else
      # warming phase (not enough runs yet)
      echo -e "$(date +'%F %T')\t$ID\t$K\t$runs\t$LL\tNA\tNA\twarming" | tee -a "$MASTER_LOG"
    fi
  done
  # Pick the best run by LL and produce labeled table for plotting
  best_line=$(awk 'NR>1{print $0}' "$BURN_TSV" | sort -k3,3g | tail -n1)
  BEST_LL=$(echo "$best_line" | awk '{print $3}')
  BEST_Q=$(echo "$best_line" | awk '{print $5}')
  BEST_F=$(echo "$best_line" | awk '{print $6}')
  echo "[K=$K] Best LL=$BEST_LL"

  cp -f "$BEST_Q" "${BEST_DIR}/${ID}.K${K}.best.qopt"
  [[ -s "$BEST_F" ]] && cp -f "$BEST_F" "${BEST_DIR}/${ID}.K${K}.best.fopt.gz" || true

  # sample IDs + Q -> labeled table
  paste "$SAMPLES_TXT" "$BEST_Q" > "${BEST_DIR}/${ID}.K${K}.best.labeled.tsv"
  echo "[K=$K] Labeled table -> ${BEST_DIR}/${ID}.K${K}.best.labeled.tsv"
done


echo "$(date) [DONE]"
