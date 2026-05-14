#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N tahoe_90th
#$ -m e
#$ -j y
#$ -l gpus=1
#$ -l gpu_c=8.0
#$ -P montilab-p
#$ -pe omp 16
#$ -l mem_per_core=16G
#$ -t 1

cd /restricted/projectnb/brcameta/projects/sig_recon/

mamba activate stack

# ------------------------------------------------
# Paths
# ------------------------------------------------

DATA_DIR="/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/tahoe/stack_subsets"
SPLIT_DIR="${DATA_DIR}/90th_splits"
CTRL_DIR="${DATA_DIR}/ctrl_comparison"

SPLITS_CSV="/restricted/projectnb/brcameta/projects/sig_recon/data/sigs/tahoe/drug_splits.csv"

OUTPUT_BASE="/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/tahoe/stack_pred/90th_perturb"

# ------------------------------------------------
# Select split file
# ------------------------------------------------

SPLIT_ARRAY=($(ls ${SPLIT_DIR}/NCI-H23_9_10th_split_*.h5ad | sort -V))
SPLIT_FILE=${SPLIT_ARRAY[$((SGE_TASK_ID-1))]}

SPLIT_NAME=$(basename ${SPLIT_FILE} .h5ad)
SPLIT_NUM=$(echo ${SPLIT_NAME} | grep -oP '\d+$')
SPLIT_COL="split_${SPLIT_NUM}"

echo "Processing ${SPLIT_NAME}"
echo "Using split column ${SPLIT_COL}"

OUTPUT_DIR="${OUTPUT_BASE}/${SPLIT_NAME}"
mkdir -p "${OUTPUT_DIR}"

# ------------------------------------------------
# Build drug list from split file
# ------------------------------------------------

mapfile -t DRUGS_ARRAY < <(awk -F',' -v col="${SPLIT_COL}" '
NR==1 {
  for(i=1;i<=NF;i++){
    gsub(/"/,"",$i)
    if($i==col) col_idx=i
  }
  next
}
{
  gsub(/"/,"",$col_idx)
  if(tolower($col_idx) ~ /^(true|True|TRUE)$/){
    gsub(/"/,"",$1)
    print $1
  }
}
' "${SPLITS_CSV}")

# always include control
DRUGS_ARRAY+=("DMSO_TF")

# ------------------------------------------------
# Build normalized set of already generated drugs
# ------------------------------------------------

normalize() {
  echo "$1" \
    | tr '[:upper:]' '[:lower:]' \
    | sed -E 's/[[:space:]]+/_/g' \
    | sed -E 's/[^a-z0-9_]+/_/g' \
    | sed -E 's/_+/_/g' \
    | sed -E 's/^_|_$//g'
}

declare -A DONE_MAP

for f in "${OUTPUT_DIR}"/*; do
  [ -e "$f" ] || continue
  base=$(basename "$f")
  base="${base%.*}"
  norm=$(normalize "$base")
  DONE_MAP["$norm"]=1
done

# ------------------------------------------------
# Remove already generated outputs
# ------------------------------------------------

ORIG_COUNT=${#DRUGS_ARRAY[@]}

FILTERED_DRUGS=()

for drug in "${DRUGS_ARRAY[@]}"; do
  norm_drug=$(normalize "$drug")

  if [[ -z "${DONE_MAP[$norm_drug]}" ]]; then
    FILTERED_DRUGS+=("$drug")
  else
    echo "Skipping $drug (already generated)"
  fi
done

DRUGS_ARRAY=("${FILTERED_DRUGS[@]}")
FILTERED_COUNT=${#DRUGS_ARRAY[@]}

echo "Drugs already generated in OUTPUT_DIR: ${DONE_COUNT}"
echo "Total drugs from split: ${ORIG_COUNT}"
echo "Drugs remaining to generate: ${FILTERED_COUNT}"
echo "-------------------------------------"

if [ ${#DRUGS_ARRAY[@]} -eq 0 ]; then
  echo "Nothing left to process for ${SPLIT_NAME}"
  exit 0
fi

# ------------------------------------------------
# Run stack-generation
# ------------------------------------------------

stack-generation \
--checkpoint "scripts/stack/notebooks/tutorial-pred-model/bc_large_aligned.ckpt" \
--base-adata "${SPLIT_FILE}" \
--test-adata "${CTRL_DIR}/NCI-H23_other_dmso.h5ad" \
--genelist "scripts/stack/notebooks/tutorial-pred-model/basecount_1000per_15000max.pkl" \
--split-column drug \
--split-values "${DRUGS_ARRAY[@]}" \
--batch-size 16 \
--num-steps 5 \
--output-dir "${OUTPUT_DIR}" \
--prompt-ratio 0.25 \
--context-ratio 0.4

echo "Completed ${SPLIT_NAME}"