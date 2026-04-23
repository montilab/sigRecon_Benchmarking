import os
import scanpy as sc
import pandas as pd

# Paths
DATA_PATH = "/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data"
PROJ_PATH = "/restricted/projectnb/brcameta/projects/sig_recon"

# Load the full datasets
k562_all = sc.read(os.path.join(PROJ_PATH, "data/stack/sciplex/k562_all.h5ad"))
a549_all = sc.read(os.path.join(PROJ_PATH, "data/stack/sciplex/a549_all.h5ad"))
mcf7_all = sc.read(os.path.join(PROJ_PATH, "data/stack/sciplex/mcf7_all.h5ad"))

# Load drug splits
splits_df = pd.read_csv(os.path.join(PROJ_PATH, "data/sigs/sciplex/drug_splits.csv"))

# Output directory
output_dir = os.path.join(PROJ_PATH, "data/stack/sciplex/splits")
os.makedirs(output_dir, exist_ok=True)

# Get perturbed (non-Vehicle) cells for each cell line
def get_perturbed(adata, drug_col='product_name'):
    """Get cells where drug_col != 'Vehicle'"""
    return adata[adata.obs[drug_col] != 'Vehicle'].copy()

# Get available drugs and perturbed cells
k562_pert = get_perturbed(k562_all)
a549_pert = get_perturbed(a549_all)
mcf7_pert = get_perturbed(mcf7_all)

k562_drugs = set(k562_pert.obs['product_name'])
a549_drugs = set(a549_pert.obs['product_name'])
mcf7_drugs = set(mcf7_pert.obs['product_name'])

def create_triple_split(base_all, base_name, 
                        add1_pert, add1_name, add1_drugs,
                        add2_pert, add2_name, add2_drugs,
                        split_num, split_type, output_dir):
    """
    Create: base_all + add1_filtered + add2_filtered (both 10th or both 90th)
    """
    # Filter add1 and add2 to specified drugs
    add1_filtered = add1_pert[add1_pert.obs['product_name'].isin(add1_drugs)].copy()
    add2_filtered = add2_pert[add2_pert.obs['product_name'].isin(add2_drugs)].copy()
    
    if add1_filtered.n_obs == 0:
        print(f"  Warning: {add1_name} has 0 cells")
    if add2_filtered.n_obs == 0:
        print(f"  Warning: {add2_name} has 0 cells")
    
    # Concatenate all three
    combined = sc.concat([base_all, add1_filtered, add2_filtered], 
                         axis=0,
                         label='dataset',
                         keys=[base_name, add1_name, add2_name],
                         index_unique='-')
    
    # Save
    out_file = os.path.join(output_dir, 
                            f"{base_name}_all_split_{split_type}_{split_num}.h5ad")
    combined.write(out_file)
    print(f"  Saved: {os.path.basename(out_file)} ({combined.n_obs} cells: "
          f"{base_all.n_obs} base + {add1_filtered.n_obs} {add1_name} + {add2_filtered.n_obs} {add2_name})")
    
    del combined, add1_filtered, add2_filtered
    import gc
    gc.collect()

# Generate all splits
print("Generating triple splits...")

for split_num in range(1, 11):
    split_col = f"split_{split_num}"
    print(f"\n=== Split {split_num} ===")
    
    # Get drugs in this split (TRUE) and not in split (FALSE)
    drugs_in_split = set(splits_df[splits_df[split_col] == True]['drug'])
    drugs_not_in_split = set(splits_df[splits_df[split_col] == False]['drug'])
    
    # Get drug subsets for each cell line
    k562_in = drugs_in_split & k562_drugs
    k562_out = drugs_not_in_split & k562_drugs
    a549_in = drugs_in_split & a549_drugs
    a549_out = drugs_not_in_split & a549_drugs
    mcf7_in = drugs_in_split & mcf7_drugs
    mcf7_out = drugs_not_in_split & mcf7_drugs
    
    # --- 10th splits (TRUE drugs) ---
    print("  10th splits:")
    
    # k562_all + a549_10th + mcf7_10th
    create_triple_split(k562_all, 'k562', 
                        a549_pert, 'a549', a549_in,
                        mcf7_pert, 'mcf7', mcf7_in,
                        split_num, '10th', output_dir)
    
    # a549_all + k562_10th + mcf7_10th
    create_triple_split(a549_all, 'a549',
                        k562_pert, 'k562', k562_in,
                        mcf7_pert, 'mcf7', mcf7_in,
                        split_num, '10th', output_dir)
    
    # mcf7_all + k562_10th + a549_10th
    create_triple_split(mcf7_all, 'mcf7',
                        k562_pert, 'k562', k562_in,
                        a549_pert, 'a549', a549_in,
                        split_num, '10th', output_dir)
    
    # --- 90th splits (FALSE drugs) ---
    print("  90th splits:")
    
    # k562_all + a549_90th + mcf7_90th
    create_triple_split(k562_all, 'k562',
                        a549_pert, 'a549', a549_out,
                        mcf7_pert, 'mcf7', mcf7_out,
                        split_num, '90th', output_dir)
    
    # a549_all + k562_90th + mcf7_90th
    create_triple_split(a549_all, 'a549',
                        k562_pert, 'k562', k562_out,
                        mcf7_pert, 'mcf7', mcf7_out,
                        split_num, '90th', output_dir)
    
    # mcf7_all + k562_90th + a549_90th
    create_triple_split(mcf7_all, 'mcf7',
                        k562_pert, 'k562', k562_out,
                        a549_pert, 'a549', a549_out,
                        split_num, '90th', output_dir)

print("\n=== Summary ===")
files = [f for f in os.listdir(output_dir) if f.endswith('.h5ad')]
print(f"Total files created: {len(files)}")
for split_type in ['10th', '90th']:
    for base in ['k562', 'a549', 'mcf7']:
        prefix = f"{base}_all_split_{split_type}_"
        count = len([f for f in files if f.startswith(prefix)])
        print(f"  {prefix}*: {count} files")