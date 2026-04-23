library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
library(sigrecon)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- "/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/srivatsan_2019"

DRUG_COL     <- "product_name"
CONTROL_NAME <- "Vehicle"

drug_splits <- read.csv(file.path(PATH, "data/sigs/sciplex/drug_splits.csv"),
                        stringsAsFactors = FALSE)

# ------------------------------------------------
# Helper: Seurat -> SummarizedExperiment (counts)
# ------------------------------------------------

prep_se <- function(seurat_obj) {
  
  counts <- GetAssayData(seurat_obj, slot = "counts")
  
  se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    colData = seurat_obj@meta.data
  )
  var_genes <- sigrecon:::rank.var.eset(se, fn = var)$values[1:10000]
  
  se[var_genes, ]
}

# ------------------------------------------------
# Load Seurat objects
# ------------------------------------------------

a549_pb <- readRDS(file.path(PATH, "data/sci_plex/a549_filtered_pb.rds"))
k562_pb <- readRDS(file.path(PATH, "data/sci_plex/k562_filtered_pb.rds"))
mcf7_pb <- readRDS(file.path(PATH, "data/sci_plex/mcf7_filtered_pb.rds"))

# ------------------------------------------------
# Convert to SummarizedExperiment
# ------------------------------------------------

a549_se <- prep_se(a549_pb)
k562_se <- prep_se(k562_pb)
mcf7_se <- prep_se(mcf7_pb)

cell_lines <- list(
  a549 = a549_se,
  k562 = k562_se,
  mcf7 = mcf7_se
)

# ------------------------------------------------
# Load signatures
# ------------------------------------------------

sigs <- list(
  a549 = lapply(sciplex.a549, function(x) x$up),
  k562 = lapply(sciplex.k562, function(x) x$up),
  mcf7 = lapply(sciplex.mcf7, function(x) x$up)
)

outdir <- file.path(PATH, "data/sigs/sciplex/projectCor")

# ------------------------------------------------
# CONTROL EXPERIMENT
# ------------------------------------------------

ctrl_se <- lapply(cell_lines, function(se) {
  ctrl <- colnames(se)[se[[DRUG_COL]] == CONTROL_NAME]
  se[, ctrl]
})

pairs <- expand.grid(source = names(cell_lines),
                     target = names(cell_lines),
                     stringsAsFactors = FALSE) %>%
  filter(source != target)

for (p in seq_len(nrow(pairs))) {
  
  src <- pairs$source[p]
  tgt <- pairs$target[p]
  
  gsva  <- projectCor(ctrl_se[[src]], sigs[[tgt]], "gsva")
  eigen <- projectCor(ctrl_se[[src]], sigs[[tgt]], "eigen")
  
  saveRDS(gsva,  file.path(outdir, paste0(tgt,"_",src,"_ctrl_gsva.rds")))
  saveRDS(eigen, file.path(outdir, paste0(tgt,"_",src,"_ctrl_eigen.rds")))
}

# ------------------------------------------------
# 1/10th and 9/10th SPLIT EXPERIMENTS
# ------------------------------------------------

for (i in 1:10) {
  
  split_col <- paste0("split_", i)
  message("Processing ", split_col)
  
  drugs_10th <- drug_splits$drug[drug_splits[[split_col]]]
  drugs_90th <- drug_splits$drug[!drug_splits[[split_col]]]
  
  for (dataset_type in c("10th","90th")) {
    
    drugs <- if (dataset_type == "10th") drugs_10th else drugs_90th
    
    subset_se <- lapply(cell_lines, function(se) {
      
      keep <- colnames(se)[
        se[[DRUG_COL]] == CONTROL_NAME |
          se[[DRUG_COL]] %in% drugs
      ]
      
      se[, keep]
    })
    
    for (p in seq_len(nrow(pairs))) {
      
      src <- pairs$source[p]
      tgt <- pairs$target[p]
      
      gsva  <- projectCor(subset_se[[src]], sigs[[tgt]], "gsva")
      eigen <- projectCor(subset_se[[src]], sigs[[tgt]], "eigen")
      
      saveRDS(gsva,
              file.path(outdir,
                        paste0(tgt,"_",src,"_",dataset_type,
                               "_gsva_split_",i,".rds")))
      
      saveRDS(eigen,
              file.path(outdir,
                        paste0(tgt,"_",src,"_",dataset_type,
                               "_eigen_split_",i,".rds")))
    }
  }
}