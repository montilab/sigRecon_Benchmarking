#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N perturb_net
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/network_prop
module load R/4.4.3
Rscript --verbose 01_networks.R
