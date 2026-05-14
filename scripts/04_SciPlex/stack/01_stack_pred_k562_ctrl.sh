#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N stack_sciplex_k562
#$ -m e
#$ -j y
#$ -l gpus=1
#$ -l gpu_c=8.0
#$ -P el-studies
#$ -pe omp 16

cd /restricted/projectnb/brcameta/projects/sig_recon/
mamba activate stack

stack-generation --checkpoint "scripts/stack/notebooks/tutorial-pred-model/bc_large_aligned.ckpt" \
--base-adata "data/stack/sciplex/k562_all.h5ad" \
--test-adata "data/stack/sciplex/a549_mcf7_ctrl.h5ad" \
--genelist "scripts/stack/notebooks/tutorial-pred-model/basecount_1000per_15000max.pkl" \
--split-column product_name \
--batch-size 16 \
--num-steps 5 \
--output-dir "/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/srivatsan_2019/stack_pred/ctrl/k562" \
--prompt-ratio 0.25 \
--context-ratio 0.4
