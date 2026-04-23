import scanpy as sc
from scipy import sparse
import anndata as ad
from adpbulk import ADPBulk
from tqdm import tqdm
import gc
import os

DATA_PATH = os.path.join(os.getenv("AGED"), "CBMrepositoryData/perturbational_data/tahoe")
SAVE_PATH = os.path.join(os.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/tahoe/pseudobulk")

for i in tqdm(range(14)):
    id = str(i + 1)
    PATH = os.path.join(DATA_PATH, f"h5ad/plate{id}_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad")
    adata = ad.read_h5ad(PATH)
    adata = adata[adata.obs.pass_filter == "full",]
    adata.raw = adata

    adpb = ADPBulk(adata, ["sample", "cell_name", "drugname_drugconc"], use_raw=True)
    pseudobulk_matrix = adpb.fit_transform()
    pseudo_meta = adpb.get_meta()

    pseudo_adata = ad.AnnData(
        X=pseudobulk_matrix.values,
        obs=pseudo_meta,
        var=adata.var.loc[pseudobulk_matrix.columns].copy()
    )
    
    output_path = os.path.join(SAVE_PATH, f"plate{id}_pseudobulk.h5ad")
    pseudo_adata.write_h5ad(output_path)

    del adata
    del pseudo_adata
    gc.collect()


data_dict = {}

for i in range(14):
    data_dict.update({f"plate_{i+1}": os.path.join(SAVE_PATH, f"plate{i+1}_pseudobulk.h5ad")})

ad.experimental.concat_on_disk(
    data_dict,
    os.path.join(SAVE_PATH, f"merged_pseudobulk.h5ad"),
    label='plate',
)


# Other changes
# evaluate_dup_counts(ad_pb)
# There are no duplicated counts.
ad_pb = ad.read_h5ad(os.path.join(SAVE_PATH, "merged_pseudobulk.h5ad"))
ad_pb.X = sparse.csr_matrix(ad_pb.X)
ad_pb.obs["drugname_drugconc"] = ad_pb.obs["drugname_drugconc"].apply(lambda x: x.replace("[","").replace("]","").replace("(","").replace(")","").replace("'",""))
ad_pb.obs["drug_name"] = ad_pb.obs["drugname_drugconc"].apply(lambda x: x.split(",")[0])
ad_pb.obs["drug_conc"] = ad_pb.obs["drugname_drugconc"].apply(lambda x: x.split(",")[1] + x.split(",")[2])
# Create unique obs_names
ad_pb.obs_names = (
    ad_pb.obs['plate'].astype(str) + '_' + 
    ad_pb.obs['cell_name'].astype(str) + "_" + 
    ad_pb.obs['drug_name'].astype(str) + "_" + 
    ad_pb.obs['drug_conc'].astype(str)
)
ad_pb.obs_names_make_unique()

# Cast categories to strings for H5AD compatibility
for col in ad_pb.obs.select_dtypes(['category']).columns:
    ad_pb.obs[col] = ad_pb.obs[col].astype(str)
for col in ad_pb.var.select_dtypes(['category']).columns:
    ad_pb.var[col] = ad_pb.var[col].astype(str)
ad_pb.write_h5ad(os.path.join(SAVE_PATH, "merged_pseudobulk.h5ad"))
print(f"Final merged pseudobulk saved.")

