#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N stack_k562_10th
#$ -m e
#$ -j y
#$ -l gpus=1
#$ -l gpu_c=8.0
#$ -P lcproject
#$ -pe omp 16
#$ -t 1-10

cd /restricted/projectnb/brcameta/projects/sig_recon/
mamba activate stack

# Paths
DATA_DIR="/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/stack_splits/perturb_seq/splits"
SPLITS_CSV="/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/sigs/perturb-seq/pb_splits.csv"

# Create array of split files and get current one
SPLIT_ARRAY=($(ls ${DATA_DIR}/k562_all_split_10th_*.h5ad | sort -V))
SPLIT_FILE=${SPLIT_ARRAY[$((SGE_TASK_ID-1))]}

# Extract split number from filename (e.g., "k562_all_split_10th_5.h5ad" -> 5)
SPLIT_NAME=$(basename ${SPLIT_FILE} .h5ad)
SPLIT_NUM=$(echo ${SPLIT_NAME} | grep -oP '\d+$')  # Extract trailing number
SPLIT_COL="split_${SPLIT_NUM}"

echo "Processing ${SPLIT_NAME} (column: ${SPLIT_COL})"
echo "Base adata: ${SPLIT_FILE}"

# Extract FALSE genes for this split into array
mapfile -t GENES_ARRAY < <(awk -F'[,\t]' -v col="${SPLIT_COL}" '
    NR==1 {
        for(i=1; i<=NF; i++) {
            gsub(/"/, "", $i)
            if($i == col) col_idx=i
        }
        next
    }
    {
        gsub(/"/, "", $col_idx)
        if(tolower($col_idx) ~ /^(false|False|FALSE|1)$/) {
            gsub(/"/, "", $1)
            print $1
        }
    }
' "${SPLITS_CSV}")

echo "Found ${#GENES_ARRAY[@]} genes for ${SPLIT_COL}"

if [ ${#GENES_ARRAY[@]} -eq 0 ]; then
    echo "ERROR: No genes found for ${SPLIT_COL}"
    exit 1
fi

echo "First 5 genes: ${GENES_ARRAY[*]:0:5}"

# Run with SPLIT_FILE as base_adata
stack-generation --checkpoint "scripts/stack/notebooks/tutorial-pred-model/bc_large_aligned.ckpt" \
--base-adata "${SPLIT_FILE}" \
--test-adata "data/stack_splits/perturb_seq/rpe1_dmso.h5ad" \
--genelist "scripts/stack/notebooks/tutorial-pred-model/basecount_1000per_15000max.pkl" \
--split-column gene \
--split-values "${GENES_ARRAY[@]}" \
--batch-size 16 \
--num-steps 5 \
--output-dir "/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/replogle_2022/stack_pred/${SPLIT_NAME}" \
--prompt-ratio 0.25 \
--context-ratio 0.4

echo "Completed ${SPLIT_NAME}"