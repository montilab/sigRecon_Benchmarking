library(argparse)
library(tidyverse)
library(anndata)
library(reticulate)
library(Matrix)
library(DESeq2)
library(doParallel)
library(foreach)
library(sigrecon)

reticulate::use_condaenv("r-sceasy", required=TRUE)
registerDoParallel(15)

parser <- ArgumentParser(description="Run DESeq2 signatures from stack splits")
parser$add_argument("--input_dir", type="character", required=TRUE)
parser$add_argument("--output_dir", type="character", required=TRUE)
parser$add_argument("--source_cell_line", type="character", required=TRUE)
parser$add_argument("--experiment", type="character", default="90th")
parser$add_argument("--ctrl_name", type="character", default="Vehicle")

args <- parser$parse_args()

DATA_PATH <- args$input_dir
SAVE_PATH <- args$output_dir
source_cell_line <- args$source_cell_line
experiment <- args$experiment
ctrl_name <- args$ctrl_name

dir.create(SAVE_PATH, showWarnings=FALSE, recursive=TRUE)

get_counts <- function(adata){
  X <- adata$X
  if(inherits(X,"scipy.sparse.csr_matrix") || inherits(X,"scipy.sparse.csc_matrix")){
    mat <- reticulate::py_to_r(X)
    mat <- as(mat,"dgCMatrix")
  } else {
    mat <- Matrix(as.matrix(X), sparse=TRUE)
  }
  mat
}

split_dirs <- paste0(DATA_PATH, "/", source_cell_line,
                     "_all_split_", experiment, "_", 1:10)

# determine target cell lines from first split
example_ctrl <- Sys.glob(file.path(split_dirs[1], paste0("*",ctrl_name,"*.h5ad")))[1]
adata_tmp <- anndata::read_h5ad(example_ctrl)
target_cell_lines <- unique(adata_tmp$obs$cell_line)

for(target_cell_line in target_cell_lines){
  
  message("Processing target cell line: ", target_cell_line)
  
  split_results <- foreach(
    split_dir = split_dirs,
    .packages=c("tidyverse","anndata","Matrix","DESeq2","reticulate"),
    .combine = 'c'
  ) %dopar% {
    
    reticulate::use_condaenv("r-sceasy", required=TRUE)
    
    split_name <- basename(split_dir)
    
    adata_filepaths <- Sys.glob(file.path(split_dir,"*.h5ad"))
    
    ctrl_filepath <- adata_filepaths[str_detect(adata_filepaths, ctrl_name)]
    drug_filepaths <- adata_filepaths[!str_detect(adata_filepaths, ctrl_name)]
    
    if(length(ctrl_filepath)==0) return(NULL)
    
    ctrl_adata <- anndata::read_h5ad(ctrl_filepath)
    ctrl_counts <- get_counts(ctrl_adata)
    
    ctrl_meta <- as.data.frame(ctrl_adata$obs) %>%
      rownames_to_column("cell_id") %>%
      filter(cell_line == target_cell_line) %>%
      mutate(condition = ctrl_name)
    
    if(nrow(ctrl_meta) < 30) return(NULL)
    
    ctrl_counts <- ctrl_counts[ctrl_meta$cell_id,,drop=FALSE]
    
    drug_results <- list()
    
    for(adata_filepath in drug_filepaths){
      
      product <- basename(adata_filepath) %>% str_remove("\\.h5ad$")
      
      pb_adata <- anndata::read_h5ad(adata_filepath)
      pb_counts <- get_counts(pb_adata)
      
      pb_meta <- as.data.frame(pb_adata$obs) %>%
        rownames_to_column("cell_id") %>%
        filter(cell_line == target_cell_line) %>%
        mutate(condition = product)
      
      if(nrow(pb_meta) < 30) next
      
      pb_counts <- pb_counts[pb_meta$cell_id,,drop=FALSE]
      
      stopifnot(all.equal(colnames(ctrl_counts), colnames(pb_counts)))
      
      counts <- rbind(ctrl_counts, pb_counts)
      meta <- bind_rows(ctrl_meta, pb_meta)
      
      meta <- meta %>%
        mutate(
          sample_id = paste0(condition,"_",replicate)
        )
      
      sample_ids <- unique(meta$sample_id)
      
      pb_matrix <- Matrix(
        0,
        nrow=ncol(counts),
        ncol=length(sample_ids),
        sparse=TRUE
      )
      
      rownames(pb_matrix) <- colnames(counts)
      colnames(pb_matrix) <- sample_ids
      
      for(s in sample_ids){
        cells <- meta %>% filter(sample_id==s) %>% pull(cell_id)
        if(length(cells)==0) next
        pb_matrix[,s] <- Matrix::colSums(counts[cells,,drop=FALSE])
      }
      
      coldata <- tibble(sample=colnames(pb_matrix)) %>%
        mutate(
          replicate=str_extract(sample,"rep[0-9]+$"),
          condition=str_remove(sample,"_rep[0-9]+$")
        )
      
      libsize <- colSums(pb_matrix)
      keep_samples <- libsize>0
      
      pb_matrix <- pb_matrix[,keep_samples,drop=FALSE]
      coldata <- coldata[keep_samples,]
      
      cond_counts <- table(coldata$condition)
      if(any(cond_counts < 2)) next
      
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
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        nbinomWaldTest(dds)
      })
      
      res <- results(dds, contrast=c("condition",product,ctrl_name))
      
      res_df <- as.data.frame(res) %>%
        rownames_to_column("gene")
      res_df$product <- product
      sig <- sig_filter_fn(res_df, 
                            product,
                            alpha = 1.1,
                            limit = 100,
                            pert_col = "product",
                            log2fc_col = "log2FoldChange",
                            pval_col = "pvalue",
                            geneid_col = "gene")
      drug_results <- c(drug_results, sig)
    }
    
    split_num <- as.numeric(sub(".*_", "", basename(split_dir)))
    list_name <- paste0("split_", split_num)
    
    setNames(list(drug_results), list_name)
  }
  
  outfile <- file.path(
    SAVE_PATH,
    paste0(source_cell_line,"_",target_cell_line,"_",experiment,"_sigs.rds")
  )
  
  saveRDS(split_results, outfile)
  
  message("Saved signatures to ", outfile)
}