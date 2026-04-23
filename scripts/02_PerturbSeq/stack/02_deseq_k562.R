library(Seurat)
library(tidyverse)
library(anndata)
library(reticulate)
reticulate::use_condaenv("r-sceasy")
options(Seurat.object.assay.version = "v5")
library(doParallel)
registerDoParallel(16)

DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/replogle_2022/k562_stack_pred/")
SAVE_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/stack/perturb_seq/")

adata_filepaths <- Sys.glob(paste0(DATA_PATH, "*.h5ad"))
ctrl_filepath <- adata_filepaths[str_detect(adata_filepaths, "non-targeting")]
adata_filepaths <- adata_filepaths[!str_detect(adata_filepaths, "non-targeting")]

k562_sigs <- foreach(adata_filepath = adata_filepaths, .combine=dplyr::bind_rows) %dopar% {

  product <- basename(adata_filepath) %>% str_remove(".h5ad")

  # Loading adata
  ctrl_adata <- anndata::read_h5ad(ctrl_filepath)
  pb_adata <- anndata::read_h5ad(adata_filepath)

  # Converting to Seurat
  ctrl_seurat <- CreateSeuratObject(counts = t(ctrl_adata$X),
                                  meta.data = ctrl_adata$obs)
  ctrl_seurat <- NormalizeData(ctrl_seurat)

  pb_seurat <- CreateSeuratObject(counts = t(pb_adata$X),
                                  meta.data = pb_adata$obs)
  pb_seurat <- NormalizeData(pb_seurat)
  pb_seurat$gene <- product

  # Merging
  pb_ctrl_seurat <- merge(ctrl_seurat, pb_seurat)
  pb_ctrl_seurat <- JoinLayers(pb_ctrl_seurat)

  # Finding Cell Ids
  pb_cell_ids <- pb_ctrl_seurat[[]] %>% rownames_to_column(var = "id") %>% dplyr::filter(gene != "non-targeting") %>% pull(id)
  ctrl_cell_ids <- pb_ctrl_seurat[[]] %>% rownames_to_column(var = "id") %>% dplyr::filter(gene == "non-targeting") %>% pull(id)

  if(length(pb_cell_ids) <= 30) {
    print(paste0("Skipping ", product, ". Fewer than 30 cells."))
    NULL
  } else {
    stopifnot(pb_cell_ids %in% colnames(pb_ctrl_seurat))
    stopifnot(ctrl_cell_ids %in% colnames(pb_ctrl_seurat))
    sig_df <- Seurat::FindMarkers(pb_ctrl_seurat,
                                  ident.1 = pb_cell_ids,
                                  ident.2 = ctrl_cell_ids,
                                  slot = "data",
                                  test.use = "MAST",
                                  latent.vars = "gem_group",
                                  logfc.threshold = 0,
                                  only.pos = FALSE,
                                  min.cells.group = 3,
                                  verbose = TRUE)
    sig_df$product <- product
    sig_df
  }
}

saveRDS(k562_sigs, file.path(SAVE_PATH, "k562_sigs.rds"))
