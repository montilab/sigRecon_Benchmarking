library(argparse)
library(tidyverse)
library(anndata)
library(reticulate)
library(Matrix)
library(DESeq2)
library(doParallel)

reticulate::use_condaenv("r-sceasy")
registerDoParallel(15)

parser <- ArgumentParser(description = "Run DESeq2 analysis on perturbation data")
parser$add_argument("--input_dir", type = "character", required = TRUE)
parser$add_argument("--output_dir", type = "character", required = TRUE)
parser$add_argument("--cell_line", type = "character", required = TRUE)
parser$add_argument("--experiment", type = "character", default = "ctrl")
parser$add_argument("--ctrl_name", type = "character", required = TRUE)

parser_args <- parser$parse_args()

DATA_PATH <- parser_args$input_dir
SAVE_PATH <- parser_args$output_dir
cell_line <- parser_args$cell_line
experiment <- parser_args$experiment
ctrl_name <- parser_args$ctrl_name

dir.create(SAVE_PATH, showWarnings = FALSE, recursive = TRUE)

adata_filepaths <- Sys.glob(paste0(DATA_PATH, "*.h5ad"))
ctrl_filepath <- adata_filepaths[str_detect(adata_filepaths, ctrl_name)]
adata_filepaths <- adata_filepaths[!str_detect(adata_filepaths, ctrl_name)]

# -----------------------
# Helper
# -----------------------

get_counts <- function(adata){
  X <- adata$X
  if(inherits(X, "scipy.sparse.csr_matrix") || inherits(X, "scipy.sparse.csc_matrix")){
    mat <- reticulate::py_to_r(X)
    mat <- as(mat, "dgCMatrix")
  } else {
    mat <- Matrix(as.matrix(X), sparse = TRUE)
  }
  return(mat) # cells x genes
}

# -----------------------
# Load control once
# -----------------------

ctrl_adata <- anndata::read_h5ad(ctrl_filepath)
ctrl_counts <- get_counts(ctrl_adata)
ctrl_meta <- as.data.frame(ctrl_adata$obs) %>%
  mutate(condition = ctrl_name) %>%
  rownames_to_column("cell_id")

# -----------------------
# Main loop
# -----------------------

all_sigs_df <- foreach(
  adata_filepath = adata_filepaths,
  .combine = bind_rows,
  .packages = c("tidyverse","anndata","Matrix","DESeq2","reticulate")
) %dopar% {
  
  reticulate::use_condaenv("r-sceasy", required=TRUE)
  
  product <- basename(adata_filepath) %>% str_remove(".h5ad")
  message("Processing ", product)
  
  pb_adata <- anndata::read_h5ad(adata_filepath)
  pb_counts <- get_counts(pb_adata)
  
  pb_meta <- as.data.frame(pb_adata$obs) %>%
    mutate(condition = product) %>%
    rownames_to_column("cell_id")
  
  if(nrow(pb_meta) <= 30){
    message("Skipping ", product, " (<30 cells)")
    return(NULL)
  }
  
  # ensure same genes
  stopifnot(all.equal(colnames(ctrl_counts), colnames(pb_counts)))
  
  if(is.null(rownames(ctrl_counts))) rownames(ctrl_counts) <- ctrl_meta$cell_id
  if(is.null(rownames(pb_counts))) rownames(pb_counts) <- pb_meta$cell_id
  
  counts <- rbind(ctrl_counts, pb_counts)
  meta <- bind_rows(ctrl_meta, pb_meta)
  
  sig_dfs <- data.frame()
  
  unique_cell_lines <- unique(meta$cell_line)
  
  for(target_cell_line in unique_cell_lines){
    
    meta_sub <- meta %>% filter(cell_line == target_cell_line)
    counts_sub <- counts[meta_sub$cell_id,]
    
    meta_sub <- meta_sub %>%
      mutate(
        condition = factor(condition),
        sample_id = paste0(condition,"_",replicate)
      )
    
    # pseudobulk
    sample_ids <- unique(meta_sub$sample_id)
    
    pb_matrix <- Matrix(
      0,
      nrow=ncol(counts_sub),
      ncol=length(sample_ids),
      sparse=TRUE
    )
    
    rownames(pb_matrix) <- colnames(counts_sub)
    colnames(pb_matrix) <- sample_ids
    
    for(s in sample_ids){
      cells <- meta_sub %>% filter(sample_id==s) %>% pull(cell_id)
      if(length(cells)==0) next
      pb_matrix[,s] <- Matrix::colSums(counts_sub[cells,,drop=FALSE])
    }
    
    coldata <- tibble(sample = colnames(pb_matrix)) %>%
      mutate(
        replicate = str_extract(sample, "rep[0-9]+$"),
        condition = str_remove(sample, "_rep[0-9]+$")
      )
    
    libsize <- colSums(pb_matrix)
    keep_samples <- libsize>0
    
    pb_matrix <- pb_matrix[,keep_samples,drop=FALSE]
    coldata <- coldata[keep_samples,]
    
    cond_counts <- table(coldata$condition)
    
    if(any(cond_counts < 2)){
      message("Skipping ",product," ",target_cell_line," (<2 replicates)")
      next
    }
    
    coldata <- as.data.frame(coldata)
    rownames(coldata) <- coldata$sample
    
    dds <- DESeqDataSetFromMatrix(
      countData = round(as.matrix(pb_matrix)),
      colData = coldata,
      design = ~ replicate + condition
    )
    
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    if(nrow(dds)==0) next
    
    dds <- tryCatch({
      DESeq(dds, quiet=TRUE)
    }, error=function(e){
      message("Falling back to gene-wise dispersions")
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersionsGeneEst(dds)
      dispersions(dds) <- mcols(dds)$dispGeneEst
      nbinomWaldTest(dds)
    })
    
    res <- results(dds, contrast=c("condition",product,ctrl_name))
    
    res_df <- as.data.frame(res) %>%
      rownames_to_column("gene") %>%
      mutate(
        product = product,
        target_cell_line = target_cell_line
      )
    
    sig_dfs <- bind_rows(sig_dfs,res_df)
  }
  
  sig_dfs
}

saveRDS(all_sigs_df,file.path(SAVE_PATH,paste0(cell_line,"_ctrl_sigs.rds")))