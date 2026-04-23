from pathlib import Path
import numpy as np
import dask.distributed as dd
import scanpy as sc
import anndata as ad
import h5py
import dask
import time
from dask.distributed import Client, LocalCluster
import os
import pybiomart

from collections import Counter
import pandas as pd
from tqdm import tqdm

sc.logging.print_header()

DATA_PATH = os.path.join(os.getenv("AGED"), "CBMrepositoryData/perturbational_data/tahoe")
SAVE_PATH = os.path.join(os.path.join(os.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/tahoe"))

def main():
    use_gpu = False

    if use_gpu:
        import rapids_singlecell as rsc
        SPARSE_CHUNK_SIZE = 100_000
        from dask_cuda import LocalCUDACluster

        import rmm
        import cupy as cp
        from cupyx.scipy import sparse as spx

        from rmm.allocators.cupy import rmm_cupy_allocator

        def set_mem():
            rmm.reinitialize(managed_memory=True)
            cp.cuda.set_allocator(rmm_cupy_allocator)
        gpus = "0,1" # comma separated like "0,1,2" for 3 gpus
        cluster = LocalCUDACluster(CUDA_VISIBLE_DEVICES=gpus)
        dask.array.register_chunk_type(spx.csr_matrix)
    else:
        SPARSE_CHUNK_SIZE = 100_000
        cluster = LocalCluster(n_workers=10)

    client = Client(cluster)

    if use_gpu:
        client.run(set_mem)
        mod = rsc
            # Add memory pool references
        cp_mempool = cp.get_default_memory_pool()
        pinned_mempool = cp.get_default_pinned_memory_pool()
    else:
        mod = sc
    
    adatas = []
    all_highly_variable_genes = []

    for i in tqdm(range(14)):
        id = str(i + 1)
        PATH = os.path.join(DATA_PATH, f"h5ad/plate{id}_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad")

        with h5py.File(PATH, "r") as f:
            adata = ad.AnnData(
                obs=ad.io.read_elem(f["obs"]),
                var=ad.io.read_elem(f["var"]),
            )
            adata.X = ad.experimental.read_elem_as_dask(
                f["X"], chunks=(SPARSE_CHUNK_SIZE, adata.shape[1])
            )
        if use_gpu:
            rsc.get.anndata_to_GPU(adata)
        # 100m filtering. According to manuscript, full filter requires cells to have at least 700 unique molecular identifiers (UMIs), less than 20% mitochondrial reads, a UMI z-score within ±3 to remove outliers in total counts, a mitochondrial percentage z-score within ±3 to remove outliers in mitochondrial content, and at least 250 genes detected.
        pass_filter_mask = adata.obs["pass_filter"] == "full"
        adata = adata[pass_filter_mask, :].copy()

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=5000)

        highly_variable_genes = set(adata.var_names[adata.var["highly_variable"]])
        all_highly_variable_genes.append(highly_variable_genes)
        adatas.append(adata)

        # Cleanup at end of iteration
        del adata
        if use_gpu:
            cp_mempool.free_all_blocks()
            pinned_mempool.free_all_blocks()

        print(f"done with preprocessing {i}")

    # select the genes appears more than two plates
    gene_counts = Counter(gene for genes in all_highly_variable_genes for gene in genes)
    # selected_genes = {gene for gene, count in gene_counts.items() if count > 2}
    # print(len(selected_genes)) #7583
    # selected_genes = {gene for gene, count in gene_counts.items() if count >= 2}
    # print(len(selected_genes)) #10755
    selected_genes = {gene for gene, count in gene_counts.items() if count >= 1}
    print(len(selected_genes)) #15625
        
    for i in tqdm(range(14)):
        id = str(i + 1)
        adata = adatas[i]
        common_genes = [g for g in adata.var_names if g in selected_genes]
        adata = adata[:, common_genes].copy()

        output_path = os.path.join(DATA_PATH, f"h5ad/preprocessed/plate{id}_filtered_preprocessed_{'gpu' if use_gpu else 'cpu'}.h5ad")
        # Ensure directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        adata.write_h5ad(output_path)

        print(f"done with writing plate {id}")
        
    data_dict = {}
    
    for i in range(14):
        data_dict.update({f"plate_{i+1}": os.path.join(DATA_PATH, f"h5ad/preprocessed/plate{i+1}_filtered_preprocessed_{'gpu' if use_gpu else 'cpu'}.h5ad")})
    
    ad.experimental.concat_on_disk(
        data_dict,
        os.path.join(DATA_PATH, f"h5ad/tahoe_plate_merged_{"gpu" if use_gpu else "cpu"}.h5ad"),
        label='plate',
    )
    print("done merging!")

if __name__ == '__main__':
    main()
