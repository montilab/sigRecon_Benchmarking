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
# mamba activate stack
# conda activate stack
conda activate /restricted/projectnb/montilab-p/personal/lkroeh/tools/stack

# Directory containing the data files
DATA_DIR="/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/tahoe/stack_subsets"
SPLIT_DIR="${DATA_DIR}/10th_splits"
CTRL_DIR="${DATA_DIR}/ctrl_comparison"

# Create array of split files
SPLIT_ARRAY=($(ls ${SPLIT_DIR}/NCI-H23_1_10th_split_*.h5ad | sort -V))

# Get the split file for this task
SPLIT_FILE=${SPLIT_ARRAY[$((SGE_TASK_ID-1))]}
SPLIT_NAME=$(basename ${SPLIT_FILE} .h5ad)

echo "Processing ${SPLIT_NAME}"

# Run stack-generation with the current split
stack-generation --checkpoint "scripts/stack/notebooks/tutorial-pred-model/bc_large_aligned.ckpt" \
--base-adata "${SPLIT_FILE}" \
--test-adata "${CTRL_DIR}/NCI-H23_other_dmso.h5ad" \
--genelist "scripts/stack/notebooks/tutorial-pred-model/basecount_1000per_15000max.pkl" \
--split-column drug \
--batch-size 16 \
--num-steps 5 \
--output-dir "/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/tahoe/stack_pred/10th_perturb/${SPLIT_NAME}" \
--prompt-ratio 0.25 \
--context-ratio 0.4