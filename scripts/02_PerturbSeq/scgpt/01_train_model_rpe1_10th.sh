#!/bin/bash -l
#$ -l h_rt=18:00:00
#$ -N scgpt_rpe1_10th
#$ -m e
#$ -j y
#$ -l gpus=1
#$ -l gpu_c=8.0
#$ -l gpu_type=L40S
#$ -t 8-9
#$ -P el-studies
#$ -pe omp 8

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/scgpt
conda activate scGPT_env
wandb login
python 01_train_model.py --data_name "rpe1" --split_num $SGE_TASK_ID --split_type "10th"