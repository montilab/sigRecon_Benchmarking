library(Seurat)
library(tidyverse)
library(anndata)
library(reticulate)
library(doParallel)
library(DESeq2)
library(SummarizedExperiment)
library(Matrix)
library(sigrecon)
registerDoParallel(27)

reticulate::use_condaenv("r-sceasy")
options(Seurat.object.assay.version = "v5")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1){
  stop("Usage: Rscript 02_deseq_stack.R <cell_line>")
}
cell_line <- args[1]
message("Running: ", cell_line)

BASE_DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/replogle_2022")
DATA_DIR <- file.path(BASE_DATA_PATH, paste0(cell_line, "_stack_pred"))
SAVE_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/sigs/perturb-seq/stack")
dir.create(SAVE_PATH, showWarnings = FALSE, recursive = TRUE)

adata_filepaths <- Sys.glob(file.path(DATA_DIR, "*.h5ad"))
ctrl_filepath <- adata_filepaths[str_detect(adata_filepaths, "non-targeting")]
adata_filepaths <- adata_filepaths[!str_detect(adata_filepaths, "non-targeting")]
if(length(ctrl_filepath) == 0) stop("Control file not found in ", DATA_DIR)

# Load control once (we'll need its counts and obs)
ctrl_adata <- anndata::read_h5ad(ctrl_filepath)

# Helper: get count matrix (genes x cells) as a sparse Matrix; handle dense or sparse X
get_counts_matrix <- function(adata){
  X <- adata$X
  if(inherits(X, "scipy.sparse.csr_matrix") || inherits(X, "scipy.sparse.csc_matrix")){
    # convert Python sparse to R Matrix
    Xr <- reticulate::py_to_r(X)
    # ensure dgCMatrix
    Xr <- as(Xr, "dgCMatrix")
  } else {
    # dense numeric matrix: transpose might be needed depending on how CreateSeuratObject was used
    Xr <- as.matrix(X)
    Xr <- Matrix(Xr, sparse = TRUE)
  }
  # anndata typically stores cells x genes in .X; we want genes x cells for DESeq2 summing by gene
  # Confirm orientation: in your earlier code you transposed when making Seurat: counts = t(adata$X)
  # So adata$X was cells x genes; after py_to_r and conversion we keep that shape and transpose downstream.
  return(Xr)
}

ctrl_counts <- get_counts_matrix(ctrl_adata)  # cells x genes
ctrl_obs <- as.data.frame(ctrl_adata$obs)

# Ensure gem_group exists
if(!"gem_group" %in% colnames(ctrl_obs)){
  stop("gem_group column not found in control obs")
}

sigs <- foreach(
  adata_filepath = adata_filepaths,
  .packages = c("reticulate","DESeq2","SummarizedExperiment","Matrix","tidyverse","stringr","tibble")
) %dopar% {
  
  reticulate::use_condaenv("r-sceasy", required = TRUE)
  product <- basename(adata_filepath) %>% str_remove(".h5ad")
  message("Processing ", product)
  
  pb_adata <- anndata::read_h5ad(adata_filepath)
  pb_counts <- get_counts_matrix(pb_adata)   # cells x genes
  pb_obs <- as.data.frame(pb_adata$obs)
  
  if(!"gem_group" %in% colnames(pb_obs)){
    stop("gem_group column not found in product obs for ", product)
  }
  
  # Build metadata with condition and gem_group for each cell
  ctrl_meta <- ctrl_obs %>% mutate(condition = "non-targeting") %>% rownames_to_column(var = "cell_id")
  pb_meta   <- pb_obs %>% mutate(condition = product) %>% rownames_to_column(var = "cell_id")
  
  # Check cells too few
  if(nrow(pb_meta) <= 30){
    message("Skipping ", product, " (<30 cells)")
    return(NULL)
  }
  
  # Combine counts and metadata (cells must have unique ids)
  if(is.null(rownames(ctrl_counts))){
    rownames(ctrl_counts) <- ctrl_meta$cell_id
  }
  if(is.null(rownames(pb_counts))){
    rownames(pb_counts) <- pb_meta$cell_id
  }
  
  # Check that genes are same in both h5ads
  stopifnot(all.equal(colnames(ctrl_counts), colnames(pb_counts)))
  
  combined_counts <- rbind(ctrl_counts, pb_counts) # cells x genes
  combined_meta <- bind_rows(ctrl_meta, pb_meta)
  
  # Now create pseudobulk by summing counts per gene within each gem_group and condition combination.
  # We'll create a sample id like paste(condition, gem_group, sep="__") and sum cells belonging to same sample.
  combined_meta <- combined_meta %>%
    mutate(gem_group = as.factor(gem_group),
           condition = as.factor(condition),
           sample_id = paste0(condition, "__", gem_group))
  
  sample_ids <- unique(combined_meta$sample_id)
  # Initialize pseudobulk matrix genes x samples
  pb_matrix <- Matrix(0, nrow = length(colnames(ctrl_counts)), ncol = length(sample_ids), sparse = TRUE)
  rownames(pb_matrix) <- colnames(ctrl_counts)
  colnames(pb_matrix) <- sample_ids
  
  # For each sample (condition x gem_group), sum counts across cells
  for(s in sample_ids){
    cell_ids <- combined_meta %>% filter(sample_id == s) %>% pull(cell_id)
    if(length(cell_ids) == 0) next
    # combined_counts is cells x genes, subset rows then colSums
    mat_sub <- combined_counts[cell_ids, , drop = FALSE]
    # ensure numeric/integer
    summed <- Matrix::colSums(mat_sub)
    pb_matrix[, s] <- summed
  }
  
  # Build colData for DESeq2
  coldata <- tibble(sample = colnames(pb_matrix)) %>%
    separate(sample, into = c("condition", "gem_group"), sep = "__", remove = FALSE) %>%
    mutate(condition = factor(condition, levels = c("non-targeting", product)),
           gem_group = factor(gem_group))
  
  # Filter out samples with zero library size
  libsize <- colSums(pb_matrix)
  keep_samples <- libsize > 0
  if(sum(keep_samples) < 2){
    message("Not enough non-zero pseudobulk samples for ", product)
    return(NULL)
  }
  pb_matrix <- pb_matrix[, keep_samples, drop = FALSE]
  coldata <- as.data.frame(coldata[keep_samples, ])
  rownames(coldata) <- coldata$sample
  
  # Require at least 3 pseudobulk replicates per condition
  cond_counts <- table(coldata$condition)
  
  if(any(cond_counts < 3)){
    message("Skipping ", product, 
            " (insufficient pseudobulk replicates: ",
            paste(names(cond_counts), cond_counts, collapse=", "), ")")
    return(NULL)
  }
  # DESeq2 dataset: counts matrix must be integer (genes x samples)
  # currently pb_matrix is genes x samples
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(pb_matrix)),
    colData = coldata,
    design = ~ gem_group + condition
  )
  
  # Pre-filter low count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  if(nrow(dds) == 0){
    message("No genes pass filter for ", product)
    return(NULL)
  }
  
  dds <- tryCatch({
    DESeq(dds, quiet = TRUE)
  }, error = function(e) {
    message("Falling back to gene-wise dispersions")
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersionsGeneEst(dds)
    dispersions(dds) <- mcols(dds)$dispGeneEst
    nbinomWaldTest(dds)
  })
  
  res <- results(dds, contrast = c("condition", product, "non-targeting"))
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column(var = "gene") %>%
    mutate(product = product)
  
  res_df
}

sigs_df <- bind_rows(sigs)
df_outfile <- paste0(cell_line, "_ctrl_sigs_df.rds")
saveRDS(sigs_df, file.path(SAVE_PATH, df_outfile))

message("Saved: ", df_outfile)