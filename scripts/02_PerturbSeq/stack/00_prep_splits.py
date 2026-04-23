import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pandas as pd
import gc

DATA_PATH = os.path.join(os.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data")
STACK_PATH = os.path.join(os.getenv("MLAB"), "projects/brcameta/projects/sig_recon/scripts/stack")

splits_df = pd.read_csv(os.path.join(DATA_PATH, "sigs/perturb-seq/pb_splits.csv"))

k562_full_data = sc.read(os.path.join(DATA_PATH, "stack_splits/perturb_seq/k562_all.h5ad"))
rpe1_full_data = sc.read(os.path.join(DATA_PATH, "stack_splits/perturb_seq/rpe1_all.h5ad"))

output_dir = os.path.join(DATA_PATH, "stack_splits/perturb_seq/splits")
os.makedirs(output_dir, exist_ok=True)

pert_col = 'gene'

# Loop through each split (1-10)
for split_num in range(2, 11):
    split_col = f"split_{split_num}"
    
    # Get genes for this split
    genes_in_split = set(splits_df[splits_df[split_col] == True]['pb'].values)
    genes_not_in_split = set(splits_df[splits_df[split_col] == False]['pb'].values)
    
    # --- VERSION 1: K562_all + RPE1 perturbed profiles where pb = TRUE (in split) ---
    
    # Filter rpe1/k562 to only perturbed profiles in the split
    k562_in_split = k562_full_data[k562_full_data.obs[pert_col].isin(genes_in_split)].copy()
    rpe1_in_split = rpe1_full_data[rpe1_full_data.obs[pert_col].isin(genes_in_split)].copy()
    
    # Concatenate
    k562_rpe1_in_split = sc.concat([k562_full_data, rpe1_in_split], 
                                    axis=0, 
                                    label='dataset', 
                                    keys=['k562', 'rpe1'],
                                    index_unique='-')
    k562_rpe1_in_split.write(os.path.join(output_dir, f"k562_all_split_10th_{split_num}.h5ad"))
    del k562_rpe1_in_split
    gc.collect()
    
    rpe1_k562_in_split = sc.concat([rpe1_full_data, k562_in_split], 
                                    axis=0, 
                                    label='dataset', 
                                    keys=['rpe1', 'k562'],
                                    index_unique='-')
    rpe1_k562_in_split.write(os.path.join(output_dir, f"rpe1_all_split_10th_{split_num}.h5ad"))
    del rpe1_k562_in_split
    gc.collect()
    
    # Save
    
    # --- VERSION 2: K562_all + RPE1 perturbed profiles where pb = FALSE (not in split) ---
    
    # Filter rpe1/k562 to only perturbed profiles in the split
    k562_in_split = k562_full_data[k562_full_data.obs[pert_col].isin(genes_not_in_split)].copy()
    rpe1_in_split = rpe1_full_data[rpe1_full_data.obs[pert_col].isin(genes_not_in_split)].copy()
    
    # Concatenate
    k562_rpe1_in_split = sc.concat([k562_full_data, rpe1_in_split], 
                                    axis=0, 
                                    label='dataset', 
                                    keys=['k562', 'rpe1'],
                                    index_unique='-')
    k562_rpe1_in_split.write(os.path.join(output_dir, f"k562_all_split_90th_{split_num}.h5ad"))
    del k562_rpe1_in_split
    gc.collect()
    
    rpe1_k562_in_split = sc.concat([rpe1_full_data, k562_in_split], 
                                    axis=0, 
                                    label='dataset', 
                                    keys=['rpe1', 'k562'],
                                    index_unique='-')
    
    rpe1_k562_in_split.write(os.path.join(output_dir, f"rpe1_all_split_90th_{split_num}.h5ad"))
    del rpe1_k562_in_split
    gc.collect()
    
    # Print summary
    print(f"  Split {split_num}: {len(genes_in_split)} genes in split, {len(genes_not_in_split)} genes out")

print("\nAll splits generated successfully!")