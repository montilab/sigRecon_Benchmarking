library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
library(sigrecon)
library(doParallel)

registerDoParallel(cores = 10)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")

# ------------------------------------------------
# Signatures
# ------------------------------------------------

source_sig <- lapply(tahoe.nci_h23, function(x) x$up)

tahoe_sigs <- list(
  "A498" = tahoe.a498,
  "HCT15" = tahoe.hct15,
  "HEC-1-A" = tahoe.hec_1_a,
  "LoVo" = tahoe.lovo,
  "MIA PaCa-2" = tahoe.miapaca_2,
  "Panc 03.27" = tahoe.panc03.27,
  "SNU-1" = tahoe.snu_1,
  "SNU-423" = tahoe.snu_423,
  "SW48" = tahoe.sw48
)

# ------------------------------------------------
# Load pseudobulk
# ------------------------------------------------

seurat_obj_f <- readRDS(file.path(PATH,"data/tahoe/pseudobulk/merged_pseudobulk_filtered.rds"))

celllines <- read.csv(file.path(PATH,"data/sigs/tahoe/cell_lines.csv"))$x
celllines <- celllines[2:10]
drug_splits <- read.csv(file.path(PATH,"data/sigs/tahoe/drug_splits.csv"))

drug_subset <- c(drug_splits$drug,"DMSO_TF")

seurat_obj_f_sub <- subset(
  seurat_obj_f,
  (cell_name %in% celllines) &
    (drug_name %in% drug_subset)
)

# ------------------------------------------------
# Precompute control datasets
# ------------------------------------------------

cellline_se <- list()

for(cell_line in celllines){
  
  seurat_cellline <- subset(seurat_obj_f_sub, cell_name == cell_line)
  
  ctrl_samples <- colnames(seurat_cellline)[seurat_cellline$drug_name == "DMSO_TF"]
  
  seurat_cellline <- seurat_cellline[,ctrl_samples]
  
  seurat_cellline <- FindVariableFeatures(seurat_cellline,nfeatures=10000)
  var_feats <- VariableFeatures(seurat_cellline)
  
  seurat_cellline <- NormalizeData(seurat_cellline)
  seurat_cellline <- seurat_cellline[var_feats,]
  cellline_se[[cell_line]] <- as.SingleCellExperiment(seurat_cellline)
}

# ------------------------------------------------
# CONTROL experiment
# ------------------------------------------------

foreach(target = names(tahoe_sigs)) %dopar% {
  
  target_se <- cellline_se[[target]]
  
  gsva  <- projectCor(target_se,source_sig,"gsva")
  eigen <- projectCor(target_se,source_sig,"eigen")
  
  saveRDS(gsva,
          file.path(PATH,
                    paste0("data/sigs/tahoe/projectCor/nci_h23_",target,"_ctrl_gsva.rds")))
  
  saveRDS(eigen,
          file.path(PATH,
                    paste0("data/sigs/tahoe/projectCor/nci_h23_",target,"_ctrl_eigen.rds")))
}

# ------------------------------------------------
# SPLIT experiments
# ------------------------------------------------

foreach(target = names(tahoe_sigs)) %dopar% {
  
  target_sig <- lapply(tahoe_sigs[[target]],function(x)x$up)
  
  seurat_cellline <- subset(seurat_obj_f_sub,cell_name == target)
  
  gsva_10th  <- list()
  eigen_10th <- list()
  
  gsva_90th  <- list()
  eigen_90th <- list()
  
  for(i in 1:10){
    
    split_col <- paste0("split_",i)
    
    drugs_10th <- drug_splits$drug[drug_splits[[split_col]]]
    drugs_90th <- drug_splits$drug[!drug_splits[[split_col]]]
    
    # -------- 1/10th --------
    
    samples <- colnames(seurat_cellline)[      seurat_cellline$drug_name == "DMSO_TF" |      seurat_cellline$drug_name %in% drugs_10th    ]
    
    seurat_tmp <- seurat_cellline[,samples]
    
    seurat_tmp <- FindVariableFeatures(seurat_tmp,nfeatures=10000)
    var_feats <- VariableFeatures(seurat_tmp)
    
    seurat_tmp <- NormalizeData(seurat_tmp)
    seurat_tmp <- seurat_tmp[var_feats,]
    se_tmp <- as.SingleCellExperiment(seurat_tmp)
    
    gsva_10th[[paste0("split_",i)]]  <- projectCor(se_tmp,source_sig,"gsva")
    eigen_10th[[paste0("split_",i)]] <- projectCor(se_tmp,source_sig,"eigen")
    
    # -------- 9/10th --------
    
    samples <- colnames(seurat_cellline)[
      seurat_cellline$drug_name == "DMSO_TF" |
        seurat_cellline$drug_name %in% drugs_90th
    ]
    
    seurat_tmp <- seurat_cellline[,samples]
    
    seurat_tmp <- FindVariableFeatures(seurat_tmp,nfeatures=10000)
    var_feats <- VariableFeatures(seurat_tmp)
    
    seurat_tmp <- NormalizeData(seurat_tmp)
    seurat_tmp <- seurat_tmp[var_feats,]
    se_tmp <- as.SingleCellExperiment(seurat_tmp)
    
    gsva_90th[[paste0("split_",i)]]  <- projectCor(se_tmp,source_sig,"gsva")
    eigen_90th[[paste0("split_",i)]] <- projectCor(se_tmp,source_sig,"eigen")
  }
  
  saveRDS(gsva_10th,
          file.path(PATH,
                    paste0("data/sigs/tahoe/projectCor/nci_h23_",target,"_10th_gsva.rds")))
  
  saveRDS(eigen_10th,
          file.path(PATH,
                    paste0("data/sigs/tahoe/projectCor/nci_h23_",target,"_10th_eigen.rds")))
  
  saveRDS(gsva_90th,
          file.path(PATH,
                    paste0("data/sigs/tahoe/projectCor/nci_h23_",target,"_90th_gsva.rds")))
  
  saveRDS(eigen_90th,
          file.path(PATH,
                    paste0("data/sigs/tahoe/projectCor/nci_h23_",target,"_90th_eigen.rds")))
}