library(argparse)
library(tidyverse)
library(anndata)
library(reticulate)
library(Matrix)
library(DESeq2)
library(doParallel)
library(foreach)
library(sigrecon)
library(mori)

reticulate::use_condaenv("r-sceasy", required=TRUE)
registerDoParallel(10)

parser <- ArgumentParser()
parser$add_argument("--input_dir", required=TRUE)
parser$add_argument("--output_dir", required=TRUE)
parser$add_argument("--experiment", default="10th_perturb")
parser$add_argument("--ctrl_name", default="DMSO_TF")

args <- parser$parse_args()

DATA_PATH <- args$input_dir
SAVE_PATH <- args$output_dir
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

split_dirs <- list.dirs(DATA_PATH, recursive=FALSE)

example_ctrl <- Sys.glob(file.path(split_dirs[1], paste0("*",ctrl_name,"*.h5ad")))[1]
adata_tmp <- anndata::read_h5ad(example_ctrl)

cell_lines <- unique(adata_tmp$obs$cell_name)

for(target_cell_line in cell_lines){
  
  message("Processing cell line: ", target_cell_line)
  
  split_results <- foreach(
    split_dir = split_dirs,
    .packages=c("tidyverse","anndata","Matrix","DESeq2","reticulate"),
    .combine='c'
  ) %dopar% {
    
    reticulate::use_condaenv("r-sceasy", required=TRUE)
    
    split_name <- basename(split_dir)
    
    adata_files <- Sys.glob(file.path(split_dir,"*.h5ad"))
    
    ctrl_file <- adata_files[str_detect(adata_files, ctrl_name)]
    drug_files <- adata_files[!str_detect(adata_files, ctrl_name)]
    
    if(length(ctrl_file)==0) return(NULL)
    
    ctrl_adata <- anndata::read_h5ad(ctrl_file)
    ctrl_counts <- get_counts(ctrl_adata)
    
    ctrl_meta <- as.data.frame(ctrl_adata$obs) %>%
      rownames_to_column("cell_id") %>%
      filter(cell_name == target_cell_line) %>%
      mutate(condition = ctrl_name)
    
    if(nrow(ctrl_meta) < 30) return(NULL)
    
    ctrl_counts <- share(ctrl_counts[ctrl_meta$cell_id,,drop=FALSE])
    
    drug_results <- list()
    
    for(drug_file in drug_files){
      
      drug_name <- basename(drug_file) %>% str_remove("\\.h5ad$")
      
      pb_adata <- anndata::read_h5ad(drug_file)
      pb_counts <- get_counts(pb_adata)
      
      pb_meta <- as.data.frame(pb_adata$obs) %>%
        rownames_to_column("cell_id") %>%
        filter(cell_name == target_cell_line) %>%
        mutate(condition = drug_name)
      
      if(nrow(pb_meta) < 30) next
      
      pb_counts <- pb_counts[pb_meta$cell_id,,drop=FALSE]
      
      stopifnot(all.equal(colnames(ctrl_counts), colnames(pb_counts)))
      
      counts <- rbind(ctrl_counts, pb_counts)
      meta <- bind_rows(ctrl_meta, pb_meta)
      
      meta <- meta %>%
        mutate(sample_id=paste0(sample,"_",cell_name, "_[", condition, "]_",plate))
      
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
      
      coldata <- tibble(sample = colnames(pb_matrix)) %>%
        mutate(
          cell_name = str_match(sample, "^smp_[0-9]+_([^_]+)")[,2],
          drug_name = str_match(sample, "\\[([^]]+)\\]")[,2],
          plate = str_match(sample, "plate_(\\d+)")[,2]
        )
      
      keep <- colSums(pb_matrix)>0
      pb_matrix <- pb_matrix[,keep,drop=FALSE]
      coldata <- coldata[keep,]
      
      if(any(table(coldata$drug_name) < 2)) next
      
      coldata <- as.data.frame(coldata)
      rownames(coldata) <- coldata$sample
      
      dds <- DESeqDataSetFromMatrix(
        countData=round(as.matrix(pb_matrix)),
        colData=coldata,
        design=~ plate + drug_name
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
      
      res <- results(dds, contrast=c("drug_name",drug_name,ctrl_name))
      
      res_df <- as.data.frame(res) %>%
        rownames_to_column("gene")
      
      res_df$drug <- drug_name
      
      sig <- sig_filter_fn(
        res_df,
        drug_name,
        alpha=1.1,
        limit=100,
        pert_col="drug",
        log2fc_col="log2FoldChange",
        pval_col="pvalue",
        geneid_col="gene"
      )
      
      drug_results <- c(drug_results, sig)
      
      rm(pb_adata, pb_counts, pb_meta, pb_matrix, coldata)
      gc()
    }
    split_num <- as.numeric(sub(".*_", "", split_name))
    setNames(list(drug_results), paste0("split_", split_num))
  }
  
  outfile <- file.path(
    SAVE_PATH,
    paste0("tahoe_",target_cell_line,"_",experiment,"_sigs.rds")
  )
  
  saveRDS(split_results, outfile)
  
  message("Saved ", outfile)
}