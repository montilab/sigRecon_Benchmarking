#!/bin/bash -l
#$ -l h_rt=16:00:00
#$ -N deseq_tahoe_ctrl
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 28

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/06_Tahoe/stack
module load R/4.4.3
Rscript --verbose 02_deseq_ctrl.R \
  --input_dir "/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/tahoe/stack_pred/ctrl_perturb/NCI-H23" \
  --output_dir "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/sigs/tahoe/stack/" \
  --ctrl_name "DMSO_TF"
