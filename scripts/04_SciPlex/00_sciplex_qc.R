library(tidyverse)
library(data.table)
library(cowplot)
library(Seurat)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Matrix)
options(Seurat.object.assay.version = 'v5')

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/srivatsan_2019")

# Define outlier detection function (MAD-based)
is_outlier <- function(seurat_obj, metric, nmads = 5) {
  values <- seurat_obj@meta.data[[metric]]
  median_val <- median(values, na.rm = TRUE)
  mad_val <- mad(values, na.rm = TRUE)
  
  # Handle case where MAD is 0
  if (mad_val == 0 || is.na(mad_val)) {
    return(rep(FALSE, length(values)))
  }
  
  lower <- median_val - nmads * mad_val
  upper <- median_val + nmads * mad_val
  
  outliers <- (values < lower | values > upper)
  outliers[is.na(outliers)] <- FALSE
  
  return(outliers)
}

qc_fn <- function(seurat_obj) {
  # Convert Ensembl IDs to gene symbols using org.Hs.eg.db
  ensembl_ids <- rownames(seurat_obj)
  
  # Map Ensembl to Symbol and Chromosome
  suppressWarnings({
    gene_symbols <- mapIds(org.Hs.eg.db,
                           keys = ensembl_ids,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
    
    chromosomes <- mapIds(org.Hs.eg.db,
                          keys = ensembl_ids,
                          column = "CHR",
                          keytype = "ENSEMBL",
                          multiVals = "first")
  })
  
  # Add to feature metadata
  seurat_obj[["RNA"]]@meta.features$gene_symbol <- gene_symbols[rownames(seurat_obj)]
  seurat_obj[["RNA"]]@meta.features$chromosome <- chromosomes[rownames(seurat_obj)]
  
  # Mark gene types
  seurat_obj[["RNA"]]@meta.features$mt <- grepl("^MT-", seurat_obj[["RNA"]]@meta.features$gene_symbol, ignore.case = TRUE) |
    seurat_obj[["RNA"]]@meta.features$chromosome == "MT"
  seurat_obj[["RNA"]]@meta.features$ribo <- grepl("^RPS|^RPL", seurat_obj[["RNA"]]@meta.features$gene_symbol)
  seurat_obj[["RNA"]]@meta.features$hb <- grepl("^HB[^P]", seurat_obj[["RNA"]]@meta.features$gene_symbol)
  
  # Get data matrix (keep as sparse)
  data_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  
  # Calculate total counts (sparse-aware)
  total_counts <- Matrix::colSums(data_matrix)
  seurat_obj$total_counts <- total_counts
  seurat_obj$n_genes_by_counts <- Matrix::colSums(data_matrix > 0)
  
  # Calculate MT percentage manually
  mt_genes <- rownames(seurat_obj)[which(seurat_obj[["RNA"]]@meta.features$mt & !is.na(seurat_obj[["RNA"]]@meta.features$mt))]
  if (length(mt_genes) > 0) {
    mt_counts <- Matrix::colSums(data_matrix[mt_genes, , drop = FALSE])
    seurat_obj$pct_counts_mt <- (mt_counts / total_counts) * 100
  } else {
    warning("No mitochondrial genes found")
    seurat_obj$pct_counts_mt <- 0
  }
  
  # Ribosomal percentage
  ribo_genes <- rownames(seurat_obj)[which(seurat_obj[["RNA"]]@meta.features$ribo & !is.na(seurat_obj[["RNA"]]@meta.features$ribo))]
  if (length(ribo_genes) > 0) {
    ribo_counts <- Matrix::colSums(data_matrix[ribo_genes, , drop = FALSE])
    seurat_obj$pct_counts_ribo <- (ribo_counts / total_counts) * 100
  } else {
    seurat_obj$pct_counts_ribo <- 0
  }
  
  # Hemoglobin genes
  hb_genes <- rownames(seurat_obj)[which(seurat_obj[["RNA"]]@meta.features$hb & !is.na(seurat_obj[["RNA"]]@meta.features$hb))]
  if (length(hb_genes) > 0) {
    hb_counts <- Matrix::colSums(data_matrix[hb_genes, , drop = FALSE])
    seurat_obj$pct_counts_hb <- (hb_counts / total_counts) * 100
  } else {
    seurat_obj$pct_counts_hb <- 0
  }
  
  # Log1p transforms for outlier detection
  seurat_obj$log1p_total_counts <- log1p(seurat_obj$total_counts)
  seurat_obj$log1p_n_genes_by_counts <- log1p(seurat_obj$n_genes_by_counts)
  
  # Summary statistics
  cat(sprintf("Min total counts: %.2f\n", min(seurat_obj$total_counts)))
  cat(sprintf("Max total counts: %.2f\n", max(seurat_obj$total_counts)))
  cat(sprintf("Max pct_counts_mt: %.2f\n", max(seurat_obj$pct_counts_mt)))
  
  # Identify outliers
  seurat_obj$outlier <- (
    is_outlier(seurat_obj, "log1p_total_counts", 5) |
    is_outlier(seurat_obj, "log1p_n_genes_by_counts", 5) 
  )
  
  cat("Outlier counts:\n")
  print(table(seurat_obj$outlier, useNA = "always"))
  
  seurat_obj$mt_outlier <- is_outlier(seurat_obj, "pct_counts_mt", 3)
  
  cat("MT outlier counts:\n")
  print(table(seurat_obj$mt_outlier, useNA = "always"))
  
  # Print cell counts
  cat(sprintf("Total number of cells: %d\n", ncol(seurat_obj)))
  
  # Filter cells
  seurat_obj_filtered <- subset(seurat_obj, 
                                 subset = outlier == FALSE & mt_outlier == FALSE)
  
  cat(sprintf("Number of cells after filtering: %d\n", ncol(seurat_obj_filtered)))
  return(seurat_obj_filtered)
}

# Process each dataset
mcf7 <- readRDS(file.path(DATA_PATH, "mcf7.rds"))
mcf7_filtered <- qc_fn(mcf7)
saveRDS(mcf7_filtered, file.path(PATH, "data/sci_plex/mcf7_filtered.rds"))
rm(mcf7, mcf7_filtered)
gc()

a549 <- readRDS(file.path(DATA_PATH, "a549.rds"))
a549_filtered <- qc_fn(a549)
saveRDS(a549_filtered, file.path(PATH, "data/sci_plex/a549_filtered.rds"))
rm(a549, a549_filtered)
gc()

k562 <- readRDS(file.path(DATA_PATH, "k562.rds"))
k562_filtered <- qc_fn(k562)
saveRDS(k562_filtered, file.path(PATH, "data/sci_plex/k562_filtered.rds"))
rm(k562, k562_filtered)
gc()