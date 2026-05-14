#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N stack_sci_mcf7_10
#$ -m e
#$ -j y
#$ -l gpus=1
#$ -l gpu_c=8.0
#$ -P montilab-p
#$ -pe omp 16
#$ -l mem_per_core=16G
#$ -t 1-10

cd /restricted/projectnb/brcameta/projects/sig_recon/
module load miniconda
conda activate stack

# Paths
DATA_DIR="/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/stack_splits/sciplex/splits"
SPLITS_CSV="/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/sigs/sciplex/drug_splits.csv"

# Create array of split files and get current one
SPLIT_ARRAY=($(ls ${DATA_DIR}/mcf7_all_split_10th_*.h5ad | sort -V))
SPLIT_FILE=${SPLIT_ARRAY[$((SGE_TASK_ID-1))]}

# Extract split number from filename (e.g., "mcf7_all_split_10th_5.h5ad" -> 5)
SPLIT_NAME=$(basename ${SPLIT_FILE} .h5ad)
SPLIT_NUM=$(echo ${SPLIT_NAME} | grep -oP '\d+$')  # Extract trailing number
SPLIT_COL="split_${SPLIT_NUM}"

echo "Processing ${SPLIT_NAME} (column: ${SPLIT_COL})"
echo "Base adata: ${SPLIT_FILE}"

# Output directory
OUTPUT_DIR="/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/srivatsan_2019/stack_pred/10th/${SPLIT_NAME}"
mkdir -p "${OUTPUT_DIR}"

# Extract FALSE DRUGS for this split into array
mapfile -t DRUGS_ARRAY < <(awk -F'[,\t]' -v col="${SPLIT_COL}" '
    NR==1 {
        for(i=1; i<=NF; i++) {
            gsub(/"/, "", $i)
            if($i == col) col_idx=i
        }
        next
    }
    {
        gsub(/"/, "", $col_idx)
        if(tolower($col_idx) ~ /^(false|False|FALSE)$/) {
            gsub(/"/, "", $1)
            print $1
        }
    }
' "${SPLITS_CSV}")

# Append control perturbation
DRUGS_ARRAY+=("Vehicle")

# ------------------------------------------------
# Remove DRUGS already processed
# ------------------------------------------------

FILTERED_DRUGS=()

for drug in "${DRUGS_ARRAY[@]}"; do
    if ! ls "${OUTPUT_DIR}"/*"${drug}"* >/dev/null 2>&1; then
        FILTERED_DRUGS+=("$drug")
    else
        echo "Skipping ${drug} (already exists)"
    fi
done

DRUGS_ARRAY=("${FILTERED_DRUGS[@]}")

echo "Remaining DRUGS to process: ${#DRUGS_ARRAY[@]}"

if [ ${#DRUGS_ARRAY[@]} -eq 0 ]; then
    echo "Nothing left to process for ${SPLIT_NAME}"
    exit 0
fi

echo "First 5 DRUGS: ${DRUGS_ARRAY[*]:0:5}"

# Run with SPLIT_FILE as base_adata
stack-generation --checkpoint "scripts/stack/notebooks/tutorial-pred-model/bc_large_aligned.ckpt" \
--base-adata "${SPLIT_FILE}" \
--test-adata "data/stack_splits/sciplex/a549_k562_ctrl.h5ad" \
--genelist "scripts/stack/notebooks/tutorial-pred-model/basecount_1000per_15000max.pkl" \
--split-column product_name \
--split-values "${DRUGS_ARRAY[@]}" \
--batch-size 16 \
--num-steps 5 \
--output-dir "${OUTPUT_DIR}" \
--prompt-ratio 0.25 \
--context-ratio 0.4

echo "Completed ${SPLIT_NAME}"