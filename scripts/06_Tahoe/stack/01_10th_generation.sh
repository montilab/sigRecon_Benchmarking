#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N tahoe_10th
#$ -m e
#$ -j y
#$ -l gpus=1
#$ -l gpu_c=8.0
#$ -P brcameta
#$ -pe omp 16
#$ -l mem_per_core=16G
#$ -t 1-10

cd /restricted/projectnb/brcameta/projects/sig_recon/

module load miniconda
conda activate /restricted/projectnb/montilab-p/personal/lkroeh/tools/stack

# ------------------------------------------------
# Paths
# ------------------------------------------------

DATA_DIR="/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/tahoe/stack_subsets"
SPLIT_DIR="${DATA_DIR}/10th_splits"
CTRL_DIR="${DATA_DIR}/ctrl_comparison"

SPLITS_CSV="/restricted/projectnb/brcameta/projects/sig_recon/data/sigs/tahoe/drug_splits.csv"

OUTPUT_BASE="/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/tahoe/stack_pred/10th_perturb"

# ------------------------------------------------
# Select split file
# ------------------------------------------------

SPLIT_ARRAY=($(ls ${SPLIT_DIR}/NCI-H23_1_10th_split_*.h5ad | sort -V))
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
  if(tolower($col_idx) ~ /^(false|False|FALSE)$/){
    gsub(/"/,"",$1)
    print $1
  }
}
' "${SPLITS_CSV}")

# always include control
DRUGS_ARRAY+=("DMSO_TF")

# ------------------------------------------------
# Remove already generated outputs
# ------------------------------------------------

FILTERED_DRUGS=()

for drug in "${DRUGS_ARRAY[@]}"; do
  if ! ls "${OUTPUT_DIR}"/*"${drug}"* >/dev/null 2>&1; then
    FILTERED_DRUGS+=("$drug")
  else
    echo "Skipping ${drug} (already generated)"
  fi
done

DRUGS_ARRAY=("${FILTERED_DRUGS[@]}")

echo "Remaining drugs: ${#DRUGS_ARRAY[@]}"

if [ ${#DRUGS_ARRAY[@]} -eq 0 ]; then
  echo "Nothing left to process."
  exit 0
fi

echo "First drugs: ${DRUGS_ARRAY[*]:0:5}"

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