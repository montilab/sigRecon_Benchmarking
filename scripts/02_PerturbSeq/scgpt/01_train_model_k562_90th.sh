#!/bin/bash -l
#$ -l h_rt=18:00:00
#$ -N scgpt_k562_90th
#$ -m e
#$ -j y
#$ -l gpus=1
#$ -l gpu_c=8.0
#$ -l gpu_type=L40S
#$ -t 1-5
#$ -P lcproject
#$ -pe omp 8

SPLIT_NUMS=(4 7 8 9 10)
SELECTED_SPLIT_NUM=${SPLIT_NUMS[$(($SGE_TASK_ID - 1))]}

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/scgpt
conda activate scGPT_env
wandb login
python 01_train_model.py --data_name "k562" --split_num $SELECTED_SPLIT_NUM --split_type "90th"