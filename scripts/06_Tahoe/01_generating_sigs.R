library(tidyverse)
library(Seurat)
library(anndata)
library(reticulate)
reticulate::use_condaenv("r-sceasy")
library(SeuratDisk)
library(Matrix)
library(doParallel)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
registerDoParallel(15)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
do_save <- TRUE

rm_low_rnaseq_counts <- function(
    eset,
    min_samples=NULL,
    class_id=NULL,
    assay_id=NULL,
    reads_per = 1000000,
    min_thresh=1
)
{
  ## BEGIN input checks
  stopifnot( methods::is(eset,"ExpressionSet") || methods::is(eset,"SummarizedExperiment") || methods::is(eset, "Seurat") )
  stopifnot( xor(is.null(min_samples),is.null(class_id)) )
  ## END input checks

  if ( methods::is(eset,"ExpressionSet") )
  {
    if ( !is.null(class_id)) {
      stopifnot( class_id %in% colnames(Biobase::pData(eset)) )
      groups <- Biobase::pData(eset)[,class_id]
      min_samples <- max(min_thresh,table(groups))
    }
    rpm <- colSums(Biobase::exprs(eset))/reads_per
    filter_ind <- t(apply(Biobase::exprs(eset), 1,function(x) {x>rpm}))
    filter_ind_rowsums <- apply(filter_ind, 1, sum)
    return(eset[filter_ind_rowsums >= min_samples,])
  }
  else if ( methods::is(eset,"SummarizedExperiment") )
  {
    if ( is.null(assay_id) )
      assay_id <- names(SummarizedExperiment::assays(eset))[1]
    stopifnot( assay_id %in% names(SummarizedExperiment::assays(eset)) )

    if ( !is.null(class_id)) {
      stopifnot( class_id %in% colnames(SummarizedExperiment::colData(eset)) )
      groups <- SummarizedExperiment::colData(eset)[,class_id]
      min_samples <- max(min_thresh,table(groups))
    }
    counts <- SummarizedExperiment::assays(eset)[[assay_id]]
    rpm <- colSums(counts)/reads_per
    filter_ind <- t(apply(counts, 1, function(x) {x>rpm}))
    filter_ind_rowsums <- apply(filter_ind, 1, sum)
    return(eset[filter_ind_rowsums >= min_samples,])
  }
  else if ( methods::is(eset,"Seurat") )
  {
    # Default to RNA assay if not specified
    if (is.null(assay_id)) assay_id <- "RNA"
    stopifnot(assay_id %in% names(eset@assays))
    counts <- Seurat::GetAssayData(eset, assay = assay_id, slot = "counts")
    meta <- eset@meta.data

    # If class_id is specified, use it to determine min_samples
    if (!is.null(class_id)) {
      stopifnot(class_id %in% colnames(meta))
      groups <- meta[, class_id]
      min_samples <- max(min_thresh, table(groups))
    }
    rpm <- colSums(counts) / reads_per
    filter_ind <- t(apply(counts, 1, function(x) {x > rpm}))
    filter_ind_rowsums <- rowSums(filter_ind)
    keep_genes <- which(filter_ind_rowsums >= min_samples)
    # Subset Seurat object to keep only the filtered genes
    return(subset(eset, features = rownames(counts)[keep_genes]))
  }
  else
    stop( "unrecognized oject type: ", class(eset) )
}

if (do_save) {
  # Loading data
  adata <- anndata::read_h5ad(file.path(PATH, "data/tahoe/pseudobulk/merged_pseudobulk.h5ad"))
  # Convert dgRMatrix to dgCMatrix and transpose
  mat_t <- as(adata$X, "CsparseMatrix")      # Convert to dgCMatrix
  mat_c_t <- t(mat_t)
  seurat_obj <- CreateSeuratObject(counts = mat_c_t,
                                   meta.data = adata$obs)

  saveRDS(seurat_obj, file.path(PATH, "data/tahoe/pseudobulk/merged_pseudobulk.rds"))
  
  highest_dose_filter_metadata <- seurat_obj@meta.data %>%
    tibble::rownames_to_column(var = "cell_id") %>%
    dplyr::mutate(drug_conc = as.numeric(drug_conc)) %>%
    dplyr::group_by(cell_name, drug_name) %>%
    dplyr::filter(drug_conc == max(drug_conc, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  seurat_obj_f <- seurat_obj[,highest_dose_filter_metadata$cell_id]
  seurat_obj_f <- rm_low_rnaseq_counts(seurat_obj_f, min_samples = 3)
  
  saveRDS(seurat_obj_f, file.path(PATH, "data/tahoe/pseudobulk/merged_pseudobulk_filtered.rds"))
} else {
  seurat_obj_f <- readRDS(file.path(PATH, "data/tahoe/pseudobulk/merged_pseudobulk_filtered.rds"))
}

celllines <- seurat_obj_f$cell_name %>% unique
drugs <- seurat_obj_f$drug_name %>% unique

sig_dfs <- foreach(cell = celllines, .combine = dplyr::bind_rows) %:%
  foreach (drug = drugs) %dopar% {
    print(paste(drug, cell))

    meta <- seurat_obj_f@meta.data
    meta$cell_id <- rownames(meta)
    drug_cells <- meta %>% filter(cell_name == cell & drug_name == drug)
    control_cells <- meta %>% filter(cell_name == cell & drug_name == "DMSO_TF")

    # Combine for DESeq2
    sub_cells <- c(drug_cells$cell_id, control_cells$cell_id)
    sub_meta <- meta %>% dplyr::filter(cell_id %in% sub_cells)
    sub_meta$condition <- ifelse(sub_meta$drug_name == drug, "treatment", "control")
    if(sum(sub_meta$condition == "treatment") < 2) {
      print(paste0("Skipping ", drug, " for ", cell, ". Fewer than 2 replicates."))
      return(NULL)
    } else if (sum(sub_meta$condition == "control") < 2) {
      print(paste0("Skipping ", drug, " for ", cell, ". Fewer than 2 controls."))
      return(NULL)
    }
    sub_count <- seurat_obj_f[,sub_cells]@assays$RNA$count

    # Run DESeq2
    dds <- DESeqDataSetFromMatrix(countData = sub_count, colData = sub_meta, design = ~ condition + plate)
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", "treatment", "control"))
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$cell_line <- cell
    res_df$drug <- drug
    res_df <- res_df %>%
      tidyr::drop_na() %>%
      dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
      dplyr::arrange(desc(logFC_pval))
    res_df
  }

saveRDS(sig_dfs, file.path(PATH, "data/tahoe/tahoe_deseq_dfs.rds"))

tahoe_tbls <- readRDS(file.path(PATH, "data/tahoe/tahoe_deseq_dfs.rds"))
sig_dfs <- tahoe_tbls %>% 
  dplyr::group_by(cell_line, drug) %>% 
  dplyr::group_split()

# tahoe deseq dfs has a mix of gene symbols and ensembl ids
map_ids <- function(genes) {
  is_ensembl <- grepl("^ENSG", genes)
  
  # Map only Ensembl IDs
  if (any(is_ensembl)) {
    ensembl_ids <- genes[is_ensembl]
    tryCatch(
      {
        symbols <- mapIds(org.Hs.eg.db,
                          keys = ensembl_ids,
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")
        
        # Replace NAs with original Ensembl ID
        symbols[is.na(symbols)] <- ensembl_ids[is.na(symbols)]
        
        # Insert back
        genes[is_ensembl] <- symbols
        message(paste("Mapped", length(ensembl_ids), "ensembl ids to gene symbols."))
      }, error = function(e) {
        warning("An error occurred: ", conditionMessage(e), call. = FALSE)
      }, finally = {
        return(genes)
      }
    )
  }
  return(genes)
}

tahoe_sigs <- list()
for(sig_df in sig_dfs) {
  cell_line <- unique(sig_df$cell_line) %>% as.character
  drug <- unique(sig_df$drug) %>% as.character
  result <- sig_filter_fn(sig_df, perts = drug, 
                          pert_col = "drug", 
                          log2fc_col = "log2FoldChange", 
                          pval_col = "padj", 
                          geneid_col = "gene")
  result <- lapply(result, function(x) lapply(x, function(y) map_ids(y)))
  tahoe_sigs[[cell_line]][[drug]] <- result[[drug]]
}

saveRDS(tahoe_sigs, file.path(PATH, "data/tahoe/tahoe_sigs.rds"))

pbs <- lapply(tahoe_sigs, function(x) names(x)) %>% purrr::reduce(., intersect)
tahoe_sigs_filtered <- lapply(tahoe_sigs, function(x) x[pbs])
sig_list <- list()
# Drugs need to have at least 5 up-regulated genes
for(cell_line in names(tahoe_sigs_filtered)) {
  cell_line_data <- tahoe_sigs_filtered[[cell_line]]
  cell_line_data <- cell_line_data[lapply(cell_line_data, function(x) (length(x$up) >= 5)) %>% unlist]
  sig_list[[cell_line]] <- cell_line_data
}
# Cell lines need to have at least 5 drugs
sig_list <- sig_list[lapply(sig_list, function(x) length(x) > 5) %>% unlist]
saveRDS(sig_list, file.path(PATH, "data/tahoe/tahoe_sigs_filtered.rds"))

for (cell_line in names(tahoe_sigs_filtered)) {
  cell_line_new <- str_replace_all(cell_line, "-", "_")
  cell_line_new <- str_replace_all(cell_line_new, " ", "")
  cell_line_new <- str_replace_all(cell_line_new, "/", "_")
  var_name <- paste0("tahoe.", tolower(cell_line_new))
  assign(var_name, sig_list[[cell_line]])
  eval(bquote(usethis::use_data(.(as.name(var_name)), overwrite = TRUE)))
}

