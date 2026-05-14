#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N deseq_mcf7_scip_10
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -l buyin=TRUE
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/04_SciPlex/stack
module load R/4.4.3
Rscript --verbose 02_deseq_splits.R \
  --input_dir "/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/srivatsan_2019/stack_pred/10th/" \
  --output_dir "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/sigs/sciplex/stack/" \
  --source_cell_line "mcf7" \
  --experiment "10th" \
  --ctrl_name "Vehicle"
