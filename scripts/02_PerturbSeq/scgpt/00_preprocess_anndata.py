import os
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
import anndata as ad
import pandas as pd
import pickle 
from gears import PertData

DATA_PATH = os.path.join(os.getenv("AGED"), "CBMrepositoryData/perturbational_data/replogle_2022")
SAVE_PATH = os.path.join(os.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/scgpt")

# From Perturbation Tutorial notebook
batch_size = 64
eval_batch_size = 64

# # 0. Loading raw data
k562_full_data = sc.read(os.path.join(DATA_PATH, "K562_essential_raw_singlecell_01.h5ad"))
rpe1_full_data = sc.read(os.path.join(DATA_PATH, "rpe1_raw_singlecell_01.h5ad"))
k562_full_data.raw = k562_full_data.copy()
rpe1_full_data.raw = rpe1_full_data.copy()

# 1. Normalizing count data
sc.pp.normalize_total(k562_full_data)
sc.pp.log1p(k562_full_data)
sc.pp.normalize_total(rpe1_full_data)
sc.pp.log1p(rpe1_full_data)

# 2. Creating necessary metadata columns https://github.com/snap-stanford/GEARS/blob/master/demo/data_tutorial.ipynb
k562_full_data.obs["condition"] = k562_full_data.obs["gene"].str.replace("non-targeting", "ctrl", regex=False)
k562_full_data.obs["condition"] = k562_full_data.obs["condition"].str.replace(r"^(?!ctrl$)(.*)$", r"\1+ctrl", regex=True)

rpe1_full_data.obs["condition"] = rpe1_full_data.obs["gene"].str.replace("non-targeting", "ctrl", regex=False)
rpe1_full_data.obs["condition"] = rpe1_full_data.obs["condition"].str.replace(r"^(?!ctrl$)(.*)$", r"\1+ctrl", regex=True)

k562_full_data.obs["cell_type"] = "K562"
rpe1_full_data.obs["cell_type"] = "RPE1"

# 3. Converting dense matrix to csr
k562_full_data.X = csr_matrix(k562_full_data.X, dtype=np.float32)
rpe1_full_data.X = csr_matrix(rpe1_full_data.X, dtype=np.float32)

# 4. Saving modified anndata objects to save path
k562_full_data.write_h5ad(os.path.join(SAVE_PATH, "k562/perturb_processed.h5ad"), compression="gzip")
rpe1_full_data.write_h5ad(os.path.join(SAVE_PATH, "rpe1/perturb_processed.h5ad"), compression="gzip")

# # Preparing no split (all data)
# pert_data = PertData("/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt")
# pert_data.load(data_path = "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt/rpe1")
# pert_data.prepare_split(split = "no_split", seed = i) 

# pert_data = PertData("/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt")
# pert_data.load(data_path = "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt/k562")
# pert_data.prepare_split(split = "no_split", seed = i) 

# # Preparing no test splits (90% train, 10% validation)
# for i in range(10):
#     pert_data = PertData("/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt")
#     pert_data.load(data_path = "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt/rpe1")
#     pert_data.prepare_split(split = "no_test", seed = i) 
    
#     pert_data = PertData("/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt")
#     pert_data.load(data_path = "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt/k562")
#     pert_data.prepare_split(split = "no_test", seed = i) 


# Intersection of genes only 
k562_full_data = sc.read(os.path.join(SAVE_PATH, "k562/perturb_processed.h5ad"))
rpe1_full_data = sc.read(os.path.join(SAVE_PATH, "rpe1/perturb_processed.h5ad"))

# 1. Get the feature names (gene IDs) from both AnnData objects
k562_genes = k562_full_data.var_names
rpe1_genes = rpe1_full_data.var_names

# 2. Find the intersection of these gene sets
shared_genes = list(set(k562_genes) & set(rpe1_genes))
k562_full_data_shared = k562_full_data[:, shared_genes].copy()
rpe1_full_data_shared = rpe1_full_data[:, shared_genes].copy()

k562_full_data_shared.write_h5ad(os.path.join(SAVE_PATH, "k562_shared/perturb_processed.h5ad"), compression="gzip")
rpe1_full_data_shared.write_h5ad(os.path.join(SAVE_PATH, "rpe1_shared/perturb_processed.h5ad"), compression="gzip")

k562_full_data_shared = sc.read(os.path.join(SAVE_PATH, "k562_shared/perturb_processed.h5ad"))
rpe1_full_data_shared = sc.read(os.path.join(SAVE_PATH, "rpe1_shared/perturb_processed.h5ad"))

# Preparing no split (all data)
pert_data = PertData("/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt")
# pert_data.new_data_process(dataset_name = "rpe1_shared", adata = rpe1_full_data_shared)
pert_data.load(data_path = "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt/rpe1_shared")
pert_data.prepare_split(split = "no_split", seed = 42) 

pert_data = PertData("/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt")
pert_data.new_data_process(dataset_name = "k562_shared", adata = k562_full_data_shared)
pert_data.load(data_path = "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt/k562_shared")
pert_data.prepare_split(split = "no_split", seed = 42) 

# # Preparing no test splits (90% train, 10% validation)
# for i in range(10):
#     pert_data = PertData("/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt")
#     pert_data.load(data_path = "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt/rpe1_shared")
#     pert_data.prepare_split(split = "no_test", seed = i) 
    
#     pert_data = PertData("/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt")
#     pert_data.load(data_path = "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/scgpt/k562_shared")
#     pert_data.prepare_split(split = "no_test", seed = i) 