library(argparse)
library(Seurat)
library(tidyverse)
library(anndata)
library(reticulate)
reticulate::use_condaenv("r-sceasy")
options(Seurat.object.assay.version = "v5")
library(doParallel)
registerDoParallel(16)

parser <- ArgumentParser(description = "Run DESeq2 analysis on perturbation data")
parser$add_argument("--input_dir", type = "character", required = TRUE, help = "Directory containing input data files")
parser$add_argument("--output_dir", type = "character", required = TRUE, help = "Directory to save output results")
parser$add_argument("--cell_line", type = "character", required = TRUE, help = "Cell line name")
parser$add_argument("--ctrl_name", type = "character", required = TRUE, help = "Control Perturbation Name (e.g. DMSO)")

parser_args <- parser$parse_args()
DATA_PATH <- parser_args$input_dir
SAVE_PATH <- parser_args$output_dir
cell_line <- parser_args$cell_line
ctrl_name <- parser_args$ctrl_name

adata_filepaths <- Sys.glob(paste0(DATA_PATH, "*.h5ad"))
ctrl_filepath <- adata_filepaths[str_detect(adata_filepaths, ctrl_name)]
adata_filepaths <- adata_filepaths[!str_detect(adata_filepaths, ctrl_name)]


all_sigs_df <- foreach(adata_filepath = adata_filepaths, .combine=dplyr::bind_rows) %dopar% {

  product <- basename(adata_filepath) %>% str_remove(".h5ad")

  # Loading adata
  ctrl_adata <- anndata::read_h5ad(ctrl_filepath)
  pb_adata <- anndata::read_h5ad(adata_filepath)

  # Converting to Seurat
  ctrl_seurat <- CreateSeuratObject(counts = t(ctrl_adata$X),
                                  meta.data = ctrl_adata$obs)
  ctrl_seurat <- NormalizeData(ctrl_seurat)
  ctrl_seurat$pb <- ctrl_name

  pb_seurat <- CreateSeuratObject(counts = t(pb_adata$X),
                                  meta.data = pb_adata$obs)
  pb_seurat <- NormalizeData(pb_seurat)
  pb_seurat$pb <- product 
  unique_cell_lines <- pb_seurat$cell_line %>% unique
  
  sig_dfs <- data.frame()
  for (target_cell_line in unique_cell_lines) {
    pb_subset <- subset(pb_seurat,
                        subset=(cell_line == target_cell_line))
    # Merging
    pb_ctrl_seurat <- merge(ctrl_seurat, pb_seurat)
    pb_ctrl_seurat <- JoinLayers(pb_ctrl_seurat)
    
    # Finding Cell Ids
    pb_cell_ids <- pb_ctrl_seurat[[]] %>% rownames_to_column(var = "id") %>% dplyr::filter(pb != ctrl_name) %>% pull(id)
    ctrl_cell_ids <- pb_ctrl_seurat[[]] %>% rownames_to_column(var = "id") %>% dplyr::filter(pb == ctrl_name) %>% pull(id)
    
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
                                    logfc.threshold = 0,
                                    only.pos = FALSE,
                                    min.cells.group = 3,
                                    verbose = TRUE)
      sig_df$product <- product
      sig_df$target_cell_line <- target_cell_line
      sig_dfs <- rbind(sig_dfs, sig_df)
    }
  }
  sig_dfs
}

saveRDS(all_sigs_df, file.path(SAVE_PATH, paste0(cell_line, "_sigs.rds")))
