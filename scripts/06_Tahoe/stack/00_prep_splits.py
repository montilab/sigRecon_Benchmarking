import anndata as ad
import os
import h5py
import ast
import dask.array as da
from dask.distributed import Client, LocalCluster
import numpy as np
import pandas as pd
import re

def normalize_drug_name(drug_name):
    """
    Drug names have slight differences in pseudobulk and single cell datasets.
    """
    drug_name = str(drug_name).strip()
    
    # Specific cases
    manual_mappings = {
        'Dapagliflozin ((2S)-1,2-propanediol, hydrate)': 'Dapagliflozin 2S-1'
    }
    
    if drug_name in manual_mappings:
        drug_name = manual_mappings[drug_name]
        
    # Remove parentheses and their contents if it's a salt/form
    drug_name = re.sub(r'\(([^)]+)\)', r'\1', drug_name)
        
    return drug_name

def convert_ensembl_to_gene_names(var_df, species='human'):
    """
    Convert Ensembl IDs to gene names where possible.
    
    Parameters:
    -----------
    var_df : pd.DataFrame
        The var dataframe with gene identifiers in the index
    species : str
        Species for gene name lookup ('human' or 'mouse')
    
    Returns:
    --------
    pd.DataFrame with updated index (gene names)
    """
    import mygene
    
    mg = mygene.MyGeneInfo()
    
    # Get the current gene IDs (from index)
    gene_ids = var_df.index.tolist()
    
    # Identify which are Ensembl IDs (start with ENSG for human, ENSMUSG for mouse)
    ensembl_pattern = 'ENSG' if species == 'human' else 'ENSMUSG'
    ensembl_mask = [str(g).startswith(ensembl_pattern) for g in gene_ids]
    ensembl_ids = [gene_ids[i] for i, mask in enumerate(ensembl_mask) if mask]
    
    print(f"Found {len(ensembl_ids)} Ensembl IDs out of {len(gene_ids)} total genes")
    
    if len(ensembl_ids) == 0:
        print("No Ensembl IDs found, returning original var_df")
        return var_df
    
    # Query mygene for gene symbols
    print("Querying MyGene.info for gene symbols...")
    results = mg.querymany(
        ensembl_ids,
        scopes='ensembl.gene',
        fields='symbol',
        species=species,
        returnall=True,
        verbose=False
    )
    
    # Create mapping dictionary
    ensembl_to_symbol = {}
    for result in results['out']:
        if 'symbol' in result and 'query' in result:
            ensembl_to_symbol[result['query']] = result['symbol']
        elif 'query' in result:
            # No symbol found, keep original ID
            ensembl_to_symbol[result['query']] = result['query']
    
    print(f"Successfully mapped {len(ensembl_to_symbol)} IDs")
    
    # Create new gene names list
    new_gene_names = []
    for gene_id in gene_ids:
        if str(gene_id).startswith(ensembl_pattern):
            # Use mapped symbol if available, otherwise keep original
            new_gene_names.append(ensembl_to_symbol.get(gene_id, gene_id))
        else:
            # Already a gene name, keep it
            new_gene_names.append(gene_id)
    
    # Handle duplicates by adding suffix
    gene_name_counts = pd.Series(new_gene_names).value_counts()
    duplicates = gene_name_counts[gene_name_counts > 1].index.tolist()
    
    if duplicates:
        print(f"Warning: Found {len(duplicates)} duplicate gene names, adding suffixes...")
        seen = {}
        final_names = []
        for name in new_gene_names:
            if name in duplicates:
                count = seen.get(name, 0)
                seen[name] = count + 1
                final_names.append(f"{name}_{count}" if count > 0 else name)
            else:
                final_names.append(name)
        new_gene_names = final_names
    
    # Update var_df
    var_df = var_df.copy()
    var_df.index = new_gene_names
    var_df.index.name = 'gene_name'
    
    # Store original IDs in a column for reference
    var_df['original_id'] = gene_ids
    
    return var_df


def main():
    # Paths
    PATH = os.path.join(os.getenv("MLAB"), "projects/brcameta/projects/sig_recon/")
    DATA_PATH = os.path.join(os.getenv("AGED"), "CBMrepositoryData/perturbational_data/tahoe/h5ad")
    OUTPUT_DIR = os.path.join(os.getenv("AGED"), "CBMrepositoryData/perturbational_data/tahoe/stack_subsets")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Create subdirectories for each data type
    OUTPUT_DIR_CTRL = os.path.join(OUTPUT_DIR, "ctrl_comparison")
    OUTPUT_DIR_10th = os.path.join(OUTPUT_DIR, "10th_splits")
    OUTPUT_DIR_90th = os.path.join(OUTPUT_DIR, "90th_splits")
    os.makedirs(OUTPUT_DIR_CTRL, exist_ok=True)
    os.makedirs(OUTPUT_DIR_10th, exist_ok=True)
    os.makedirs(OUTPUT_DIR_90th, exist_ok=True)

    # Optimized for 16 cores, 256 GB RAM
    SPARSE_CHUNK_SIZE = 100_000
    
    cluster = LocalCluster(
        n_workers=8,
        threads_per_worker=2,
        memory_limit='28GB',
        processes=True
    )

    # Load data
    print("Loading data...")
    with h5py.File(os.path.join(DATA_PATH, "tahoe_plate_merged_cpu.h5ad"), "r") as f:
        adata = ad.AnnData(
            obs=ad.io.read_elem(f["obs"]),
            var=ad.io.read_elem(f["var"]),
        )
        adata.X = ad.experimental.read_elem_as_dask(
            f["X"], chunks=(SPARSE_CHUNK_SIZE, adata.shape[1])
        )
    
    # Convert gene IDs to gene names
    print("\nConverting Ensembl IDs to gene names...")
    adata.var = convert_ensembl_to_gene_names(adata.var, species='human')
    print(f"First 10 gene names: {adata.var.index[:10].tolist()}")
    
    # Separating drugconc column
    adata.obs['drugconc'] = adata.obs['drugname_drugconc'].apply(
        lambda x: ast.literal_eval(x)[0][1]
    )
    
    # Load cell lines and drug splits
    cell_names = pd.read_csv(os.path.join(PATH, "data/sigs/tahoe/cell_lines.csv"))
    cell_names = cell_names.x.to_list()
    print(f"\nCell lines selected: {cell_names}")
    
    drug_splits = pd.read_csv(os.path.join(PATH, "data/sigs/tahoe/drug_splits.csv"))
    n_splits = sum(['split_' in col for col in drug_splits.columns])
    print(f"\nLoaded drug splits with {n_splits} splits")
    print(f"Total drugs in splits: {len(drug_splits)}")
    
    print("\nNormalizing drug names...")
    
    # Normalize drug names in adata and drug_splits
    ## Keep original names
    adata.obs['drug_original'] = adata.obs['drug'].copy()
    drug_splits['drug_original'] = drug_splits['drug'].copy()
    ## Apply changes
    drug_splits['drug'] = drug_splits['drug'].apply(normalize_drug_name)
    adata.obs['drug'] = adata.obs['drug'].apply(normalize_drug_name)
    
    # Verify match
    shared_drugs = set(drug_splits['drug'].values)
    adata_drugs = set(adata.obs['drug'].unique())
    
    all_present = shared_drugs.issubset(adata_drugs)
    missing = shared_drugs - adata_drugs
    
    print(f"All shared drugs present: {all_present}")
    print(f"Missing drugs: {len(missing)}")
    
    # Get all drugs from drug_splits (shared drugs)        
    shared_drugs = drug_splits['drug'].tolist()
    # ========================================================================
    # CTRL: Each cell line max conc + DMSO from all other cell lines
    # ========================================================================
#     print("\n" + "="*80)
#     print("TYPE 1: Control Comparison (Each cell line vs all others)")
#     print("="*80)

#     for cell_name in cell_names:
#         print(f"\nProcessing Type 1 for {cell_name}...")   
#         # === SUBSET 1: Max concentration + DMSO for this cell line ===
#         cell_name_mask = adata.obs['cell_name'] == cell_name
#         cell_name_obs = adata.obs[cell_name_mask]
        
#         # Find max concentration for each drug
#         max_conc_per_drug = cell_name_obs.groupby('drug', observed=True)['drugconc'].max()
        
#         # Create boolean mask for max concentration OR DMSO (only for shared drugs)
#         high_conc_mask = np.zeros(len(adata.obs), dtype=bool)
        
#         for drug, max_conc in max_conc_per_drug.items():
#             # Only include drugs that are in shared_drugs or DMSO
#             if drug not in shared_drugs and drug != 'DMSO_TF':
#                 continue
                
#             if drug == 'DMSO_TF':
#                 mask = (adata.obs['cell_name'] == cell_name) & (adata.obs['drug'] == drug)
#             else:
#                 mask = (adata.obs['cell_name'] == cell_name) & \
#                        (adata.obs['drug'] == drug) & \
#                        (adata.obs['drugconc'] == max_conc)
#             high_conc_mask |= mask.values
        
#         # Convert boolean mask to Dask array and use it for filtering
#         high_conc_mask_da = da.from_array(high_conc_mask, chunks=SPARSE_CHUNK_SIZE)
        
#         # Filter using boolean indexing
#         X_filtered = adata.X[high_conc_mask_da, :]
#         obs_subset = adata.obs[high_conc_mask].copy()
#         var_subset = adata.var.copy()
        
#         # Now compute
#         print(f"  Computing subset 1 for {cell_name}...")
#         X_computed = X_filtered.compute()
        
#         adata_subset1 = ad.AnnData(X=X_computed, obs=obs_subset, var=var_subset)
#         n_unique_drugs = len(adata_subset1.obs.drug.unique())
#         print(f" #N unique drugs {n_unique_drugs}")
#         output_path1 = os.path.join(OUTPUT_DIR_CTRL, f"{cell_name}_maxconc_dmso.h5ad")
#         adata_subset1.write_h5ad(output_path1)
#         print(f"  Written: {output_path1} ({adata_subset1.shape[0]} cells)")
        
#         del X_computed, X_filtered, obs_subset, adata_subset1
        
#         # === SUBSET 2: DMSO from all OTHER cell lines ===
#         other_cell_names_mask = (adata.obs['cell_name'] != cell_name) & \
#                                 (adata.obs['cell_name'].isin(cell_names)) & \
#                                 (adata.obs['drug'] == 'DMSO_TF')
        
#         other_cell_names_mask = other_cell_names_mask.values
#         other_mask_da = da.from_array(other_cell_names_mask, chunks=SPARSE_CHUNK_SIZE)
        
#         # Filter using boolean indexing
#         X_filtered2 = adata.X[other_mask_da, :]
#         obs_subset2 = adata.obs[other_cell_names_mask].copy()
#         var_subset2 = adata.var.copy()
        
#         # Now compute
#         print(f"  Computing subset 2 for {cell_name}...")
#         X_computed2 = X_filtered2.compute()
        
#         adata_subset2 = ad.AnnData(X=X_computed2, obs=obs_subset2, var=var_subset2)
#         n_unique_drugs = len(adata_subset2.obs.drug.unique())
#         print(f" #N unique drugs {n_unique_drugs}")
#         output_path2 = os.path.join(OUTPUT_DIR_CTRL, f"{cell_name}_other_dmso.h5ad")
#         adata_subset2.write_h5ad(output_path2)
#         print(f"  Written: {output_path2} ({adata_subset2.shape[0]} cells)")
        
#         del X_computed2, X_filtered2, obs_subset2, adata_subset2
    
    # ========================================================================
    # TYPE 2: 1/10th splits (one cell line + 1/10 from others)
    # ========================================================================
    print("\n" + "="*80)
    print("TYPE 2: 1/10th Splits")
    print("="*80)
    
    # Choose first cell line for types 2 and 3
    target_cell_line = cell_names[0]
    other_cell_lines = [c for c in cell_names if c != target_cell_line]
    
    print(f"\nTarget cell line: {target_cell_line}")
    print(f"Other cell lines: {other_cell_lines}")
    
    for split_idx in range(1, n_splits + 1):
        print(f"\nProcessing Type 2 - Split {split_idx}/{n_splits}...")
        
        # Get drugs for this split
        split_col = f'split_{split_idx}'
        drugs_in_split = drug_splits[drug_splits[split_col]]['drug'].tolist()
        print(f"  Drugs in split {split_idx}: {len(drugs_in_split)}")
        
        # Create mask for cells to include
        combined_mask = np.zeros(len(adata.obs), dtype=bool)
        
        # 1. Target cell line: max conc for all drugs in source cell line + DMSO
        target_mask = adata.obs['cell_name'] == target_cell_line
        target_obs = adata.obs[target_mask]

        max_conc_target = target_obs.groupby('drug', observed=True)['drugconc'].max()

        for drug, max_conc in max_conc_target.items():
            if drug not in shared_drugs and drug != 'DMSO_TF':
                continue
            if drug == 'DMSO_TF':
                mask = (adata.obs['cell_name'] == target_cell_line) & \
                       (adata.obs['drug'] == drug)
            else:
                mask = (adata.obs['cell_name'] == target_cell_line) & \
                       (adata.obs['drug'] == drug) & \
                       (adata.obs['drugconc'] == max_conc)

            combined_mask |= mask.values
        
        # 2. Other cell lines: max conc for drugs in this split only
        for other_cell in other_cell_lines:
            other_mask = adata.obs['cell_name'] == other_cell
            other_obs = adata.obs[other_mask]
            
            max_conc_other = other_obs.groupby('drug', observed=True)['drugconc'].max()
            
            for drug, max_conc in max_conc_other.items():
                if drug in drugs_in_split:
                    mask = (adata.obs['cell_name'] == other_cell) & \
                           (adata.obs['drug'] == drug) & \
                           (adata.obs['drugconc'] == max_conc)
                    combined_mask |= mask.values
        
        # Convert to Dask array
        combined_mask_da = da.from_array(combined_mask, chunks=SPARSE_CHUNK_SIZE)
        
        # Filter
        X_filtered = adata.X[combined_mask_da, :]
        obs_subset = adata.obs[combined_mask].copy()
        var_subset = adata.var.copy()
        
        # Compute
        print(f"  Computing split {split_idx}...")
        X_computed = X_filtered.compute()
        
        adata_split = ad.AnnData(X=X_computed, obs=obs_subset, var=var_subset)
        n_unique_drugs = len(adata_split.obs.drug.unique())
        print(f" #N unique drugs {n_unique_drugs}")
        output_path = os.path.join(OUTPUT_DIR_10th, 
                                   f"{target_cell_line}_1_10th_split_{split_idx}.h5ad")
        adata_split.write_h5ad(output_path)
        print(f"  Written: {output_path} ({adata_split.shape[0]} cells)")
        
        del X_computed, X_filtered, obs_subset, adata_split
    
    # ========================================================================
    # TYPE 3: 9/10th splits (one cell line + 9/10 from others)
    # ========================================================================
    print("\n" + "="*80)
    print("TYPE 3: 9/10th Splits")
    print("="*80)
    
    for split_idx in range(1, n_splits + 1):
        print(f"\nProcessing Type 3 - Split {split_idx}/{n_splits} (excluding split {split_idx})...")
        
        # Get drugs NOT in this split (9/10)
        split_col = f'split_{split_idx}'
        drugs_not_in_split = drug_splits[~drug_splits[split_col]]['drug'].tolist()
        print(f"  Drugs NOT in split {split_idx}: {len(drugs_not_in_split)}")
        
        # Create mask for cells to include
        combined_mask = np.zeros(len(adata.obs), dtype=bool)
        target_mask = adata.obs['cell_name'] == target_cell_line
        
        # 1. Target cell line: max conc for all drugs in source cell line + DMSO
        target_mask = adata.obs['cell_name'] == target_cell_line
        target_obs = adata.obs[target_mask]

        max_conc_target = target_obs.groupby('drug', observed=True)['drugconc'].max()

        for drug, max_conc in max_conc_target.items():
            if drug not in shared_drugs and drug != 'DMSO_TF':
                continue
            if drug == 'DMSO_TF':
                mask = (adata.obs['cell_name'] == target_cell_line) & \
                       (adata.obs['drug'] == drug)
            else:
                mask = (adata.obs['cell_name'] == target_cell_line) & \
                       (adata.obs['drug'] == drug) & \
                       (adata.obs['drugconc'] == max_conc)

            combined_mask |= mask.values
        
        # 2. Other cell lines: max conc for drugs NOT in this split
        for other_cell in other_cell_lines:
            other_mask = adata.obs['cell_name'] == other_cell
            other_obs = adata.obs[other_mask]
            
            max_conc_other = other_obs.groupby('drug', observed=True)['drugconc'].max()
            
            for drug, max_conc in max_conc_other.items():
                if drug in drugs_not_in_split:
                    mask = (adata.obs['cell_name'] == other_cell) & \
                           (adata.obs['drug'] == drug) & \
                           (adata.obs['drugconc'] == max_conc)
                    combined_mask |= mask.values
        
        # Convert to Dask array
        combined_mask_da = da.from_array(combined_mask, chunks=SPARSE_CHUNK_SIZE)
        
        # Filter
        X_filtered = adata.X[combined_mask_da, :]
        obs_subset = adata.obs[combined_mask].copy()
        var_subset = adata.var.copy()
        
        # Compute
        print(f"  Computing split {split_idx}...")
        X_computed = X_filtered.compute()
        
        adata_split = ad.AnnData(X=X_computed, obs=obs_subset, var=var_subset)
        n_unique_drugs = len(adata_split.obs.drug.unique())
        print(f" #N unique drugs {n_unique_drugs}")
        output_path = os.path.join(OUTPUT_DIR_90th, 
                                   f"{target_cell_line}_9_10th_split_{split_idx}.h5ad")
        adata_split.write_h5ad(output_path)
        print(f"  Written: {output_path} ({adata_split.shape[0]} cells)")
        
        del X_computed, X_filtered, obs_subset, adata_split

    print("\n" + "="*80)
    print("All processing complete!")
    print("="*80)
    print(f"Type 1 files: {OUTPUT_DIR_CTRL}")
    print(f"Type 2 files: {OUTPUT_DIR_10th}")
    print(f"Type 3 files: {OUTPUT_DIR_90th}")
    
    cluster.close()

if __name__ == '__main__':
    main()