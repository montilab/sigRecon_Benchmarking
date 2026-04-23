library(tidyverse)
library(Biobase)
library(GSVA)
library(sigrecon)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/drugmatrix/")

# Load signatures
drugmatrix.liver.sigs <- drugmatrix.liver
drugmatrix.kidney.sigs <- drugmatrix.kidney

kidney_sig_up <- lapply(drugmatrix.kidney.sigs, function(x) x$up)
liver_sig_up <- lapply(drugmatrix.liver.sigs, function(x) x$up)

# All data
liver_eset <- readRDS(file=file.path(DATA_PATH, "liver.rds"))
kidney_eset <- readRDS(file=file.path(DATA_PATH, "kidney.rds"))
liver_gene_symbols <- make.unique(fData(liver_eset)$`Gene Symbol`)
featureNames(liver_eset) <- liver_gene_symbols
kidney_gene_symbols <- make.unique(fData(kidney_eset)$`Gene Symbol`)
featureNames(kidney_eset) <- kidney_gene_symbols

drug_splits <- read.csv(file.path(PATH, "data/sigs/drugmatrix/drug_splits.csv"),
                        stringsAsFactors = FALSE)

# # 1. Control Experiment
# liver_control_samples <- liver_eset %>% pData %>% filter(`dose:ch1` == "0 mg/kg") %>% rownames
# kidney_control_samples <- kidney_eset %>% pData %>% filter(`dose:ch1` == "0 mg/kg") %>% rownames
# liver_eset_filtered <- liver_eset[,liver_control_samples]
# kidney_eset_filtered <- kidney_eset[,kidney_control_samples]
# 
# 
# # Filtering Features (to control only)
# liver_var_genes <- sigrecon:::rank.var.eset(liver_eset_filtered, fn = var)$values[1:10000]
# kidney_var_genes <- sigrecon:::rank.var.eset(kidney_eset_filtered, fn = var)$values[1:10000]
# ## Subsetting Esets with non-duplicated MAD genes
# liver_eset_filtered <- liver_eset_filtered[liver_var_genes,]
# kidney_eset_filtered <- kidney_eset_filtered[kidney_var_genes,]
# liver_eset_filtered <- as(liver_eset_filtered, "SummarizedExperiment")
# kidney_eset_filtered <- as(kidney_eset_filtered, "SummarizedExperiment")
# liver_eset_filtered@metadata$annotation <- NULL
# kidney_eset_filtered@metadata$annotation <- NULL
# 
# # Recontextualization
# kidney_liver_recon_gsva <- projectCor(liver_eset_filtered, kidney_sig_up, "gsva")
# kidney_liver_recon_eigen <- projectCor(liver_eset_filtered, kidney_sig_up, "eigen")
# liver_kidney_recon_gsva <- projectCor(kidney_eset_filtered, liver_sig_up, "gsva")
# liver_kidney_recon_eigen <- projectCor(kidney_eset_filtered, liver_sig_up, "eigen")
# 
# saveRDS(kidney_liver_recon_gsva, file.path(PATH, "data/sigs/drugmatrix/projectCor/kidney_liver_ctrl_gsva.rds"))
# saveRDS(kidney_liver_recon_eigen, file.path(PATH, "data/sigs/drugmatrix/projectCor/kidney_liver_ctrl_eigen.rds"))
# saveRDS(liver_kidney_recon_gsva, file.path(PATH, "data/sigs/drugmatrix/projectCor/liver_kidney_ctrl_gsva.rds"))
# saveRDS(liver_kidney_recon_eigen, file.path(PATH, "data/sigs/drugmatrix/projectCor/liver_kidney_ctrl_eigen.rds"))

# --------------------------------------------
# 2. 1/10th and 9/10th perturbation experiments
# --------------------------------------------

for (i in 1:10) {
  
  split_col <- paste0("split_", i)
  message("Processing ", split_col)
  
  # 1/10th: drugs where split_col == TRUE
  drugs_in_split <- drug_splits$drug[drug_splits[[split_col]]]
  
  # 9/10th: drugs where split_col == FALSE
  drugs_not_in_split <- drug_splits$drug[!drug_splits[[split_col]]]
  
  # ------------------------
  # 1/10th datasets
  # ------------------------
  
  liver_samples_10th <- liver_eset %>%
    pData() %>%
    filter(`dose:ch1` == "0 mg/kg" | `compound:ch1` %in% drugs_in_split) %>%
    rownames()
  
  kidney_samples_10th <- kidney_eset %>%
    pData() %>%
    filter(`dose:ch1` == "0 mg/kg" | `compound:ch1` %in% drugs_in_split) %>%
    rownames()
  
  liver_eset_10th <- liver_eset[, liver_samples_10th]
  kidney_eset_10th <- kidney_eset[, kidney_samples_10th]
  
  # feature filtering
  liver_var <- sigrecon:::rank.var.eset(liver_eset_10th, fn = var)$values[1:10000]
  kidney_var <- sigrecon:::rank.var.eset(kidney_eset_10th, fn = var)$values[1:10000]
  
  liver_eset_10th <- liver_eset_10th[liver_var,]
  kidney_eset_10th <- kidney_eset_10th[kidney_var,]
  
  liver_eset_10th <- as(liver_eset_10th, "SummarizedExperiment")
  kidney_eset_10th <- as(kidney_eset_10th, "SummarizedExperiment")
  
  liver_eset_10th@metadata$annotation <- NULL
  kidney_eset_10th@metadata$annotation <- NULL
  
  # recontextualization
  kidney_liver_10th_gsva  <- projectCor(liver_eset_10th, kidney_sig_up, "gsva")
  kidney_liver_10th_eigen <- projectCor(liver_eset_10th, kidney_sig_up, "eigen")
  
  liver_kidney_10th_gsva  <- projectCor(kidney_eset_10th, liver_sig_up, "gsva")
  liver_kidney_10th_eigen <- projectCor(kidney_eset_10th, liver_sig_up, "eigen")
  
  saveRDS(kidney_liver_10th_gsva,
          file.path(PATH, paste0("data/sigs/drugmatrix/projectCor/kidney_liver_10th_gsva_split_", i, ".rds")))
  saveRDS(kidney_liver_10th_eigen,
          file.path(PATH, paste0("data/sigs/drugmatrix/projectCor/kidney_liver_10th_eigen_split_", i, ".rds")))
  
  saveRDS(liver_kidney_10th_gsva,
          file.path(PATH, paste0("data/sigs/drugmatrix/projectCor/liver_kidney_10th_gsva_split_", i, ".rds")))
  saveRDS(liver_kidney_10th_eigen,
          file.path(PATH, paste0("data/sigs/drugmatrix/projectCor/liver_kidney_10th_eigen_split_", i, ".rds")))
  
  
  # ------------------------
  # 9/10th datasets
  # ------------------------
  
  liver_samples_90th <- liver_eset %>%
    pData() %>%
    filter(`dose:ch1` == "0 mg/kg" | `compound:ch1` %in% drugs_not_in_split) %>%
    rownames()
  
  kidney_samples_90th <- kidney_eset %>%
    pData() %>%
    filter(`dose:ch1` == "0 mg/kg" | `compound:ch1` %in% drugs_not_in_split) %>%
    rownames()
  
  liver_eset_90th <- liver_eset[, liver_samples_90th]
  kidney_eset_90th <- kidney_eset[, kidney_samples_90th]
  
  # feature filtering
  liver_var <- sigrecon:::rank.var.eset(liver_eset_90th, fn = var)$values[1:10000]
  kidney_var <- sigrecon:::rank.var.eset(kidney_eset_90th, fn = var)$values[1:10000]
  
  liver_eset_90th <- liver_eset_90th[liver_var,]
  kidney_eset_90th <- kidney_eset_90th[kidney_var,]
  
  liver_eset_90th <- as(liver_eset_90th, "SummarizedExperiment")
  kidney_eset_90th <- as(kidney_eset_90th, "SummarizedExperiment")
  
  liver_eset_90th@metadata$annotation <- NULL
  kidney_eset_90th@metadata$annotation <- NULL
  
  # recontextualization
  kidney_liver_90th_gsva  <- projectCor(liver_eset_90th, kidney_sig_up, "gsva")
  kidney_liver_90th_eigen <- projectCor(liver_eset_90th, kidney_sig_up, "eigen")
  
  liver_kidney_90th_gsva  <- projectCor(kidney_eset_90th, liver_sig_up, "gsva")
  liver_kidney_90th_eigen <- projectCor(kidney_eset_90th, liver_sig_up, "eigen")
  
  saveRDS(kidney_liver_90th_gsva,
          file.path(PATH, paste0("data/sigs/drugmatrix/projectCor/kidney_liver_90th_gsva_split_", i, ".rds")))
  saveRDS(kidney_liver_90th_eigen,
          file.path(PATH, paste0("data/sigs/drugmatrix/projectCor/kidney_liver_90th_eigen_split_", i, ".rds")))
  
  saveRDS(liver_kidney_90th_gsva,
          file.path(PATH, paste0("data/sigs/drugmatrix/projectCor/liver_kidney_90th_gsva_split_", i, ".rds")))
  saveRDS(liver_kidney_90th_eigen,
          file.path(PATH, paste0("data/sigs/drugmatrix/projectCor/liver_kidney_90th_eigen_split_", i, ".rds")))
}