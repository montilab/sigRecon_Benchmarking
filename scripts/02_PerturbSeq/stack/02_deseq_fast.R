library(doParallel)
library(foreach)
library(Seurat)
library(tidyverse)
library(anndata)
library(reticulate)

n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

reticulate::use_condaenv("r-sceasy")
options(Seurat.object.assay.version = "v5")

# ----------------------------
# Parse command line arguments
# ----------------------------

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3){
  stop("Usage: Rscript 02_deseq_stack.R <cell_line> <experiment> <split>")
}

cell_line <- args[1]
experiment <- args[2]   # 10th or 90th
split <- args[3]

message("Running: ", cell_line, " ", experiment, " split ", split)

# ----------------------------
# Paths
# ----------------------------

BASE_DATA_PATH <- file.path(
  Sys.getenv("AGED"),
  "CBMrepositoryData/perturbational_data/replogle_2022/stack_pred"
)

DATA_DIR <- file.path(
  BASE_DATA_PATH,
  paste0(cell_line, "_all_split_", experiment, "_", split)
)

SAVE_PATH <- file.path(
  Sys.getenv("MLAB"),
  "projects/brcameta/projects/sig_recon/data/stack/perturb_seq"
)

dir.create(SAVE_PATH, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Locate files
# ----------------------------

adata_filepaths <- Sys.glob(file.path(DATA_DIR, "*.h5ad"))

ctrl_filepath <- adata_filepaths[str_detect(adata_filepaths, "non-targeting")]
adata_filepaths <- adata_filepaths[!str_detect(adata_filepaths, "non-targeting")]

if(length(ctrl_filepath) == 0){
  stop("Control file not found in ", DATA_DIR)
}

# ----------------------------
# Load control once
# ----------------------------

ctrl_adata <- anndata::read_h5ad(ctrl_filepath)

ctrl_seurat <- CreateSeuratObject(
  counts = t(ctrl_adata$X),
  meta.data = ctrl_adata$obs
)

ctrl_seurat <- NormalizeData(ctrl_seurat)

# ----------------------------
# Run differential expression
# ----------------------------

sigs <- foreach(
  adata_filepath = adata_filepaths,
  .packages = c("Seurat","tidyverse","anndata","stringr","tibble")
) %dopar% {
  
  product <- basename(adata_filepath) %>% str_remove(".h5ad")
  
  message("Processing ", product)
  
  pb_adata <- anndata::read_h5ad(adata_filepath)
  
  pb_seurat <- CreateSeuratObject(
    counts = t(pb_adata$X),
    meta.data = pb_adata$obs
  )
  
  pb_seurat <- NormalizeData(pb_seurat)
  pb_seurat$gene <- product
  
  pb_ctrl_seurat <- merge(ctrl_seurat, pb_seurat)
  pb_ctrl_seurat <- JoinLayers(pb_ctrl_seurat)
  
  meta <- pb_ctrl_seurat[[]] %>% rownames_to_column("id")
  
  pb_cell_ids <- meta %>% filter(gene != "non-targeting") %>% pull(id)
  ctrl_cell_ids <- meta %>% filter(gene == "non-targeting") %>% pull(id)
  
  if(length(pb_cell_ids) <= 30){
    message("Skipping ", product, " (<30 cells)")
    return(NULL)
  }
  
  sig_df <- Seurat::FindMarkers(
    pb_ctrl_seurat,
    ident.1 = pb_cell_ids,
    ident.2 = ctrl_cell_ids,
    slot = "data",
    test.use = "MAST",
    latent.vars = "gem_group",
    logfc.threshold = 0,
    only.pos = FALSE,
    min.cells.group = 3,
    verbose = FALSE
  )
  
  sig_df$product <- product
  sig_df
}

sigs_all <- bind_rows(sigs)

# ----------------------------
# Save
# ----------------------------

outfile <- paste0(cell_line, "_", experiment, "_split_", split, "_sigs.rds")

saveRDS(
  sigs_all,
  file.path(SAVE_PATH, outfile)
)

message("Saved: ", outfile)