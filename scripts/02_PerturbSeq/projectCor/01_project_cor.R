library(tidyverse)
library(SummarizedExperiment)
library(GSVA)
library(sigrecon)
library(doParallel)
registerDoParallel(cores=10)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
# ------------------------------------------------
# Load signatures
# ------------------------------------------------

k562_sig_up <- lapply(perturbseq.k562, function(x) x$up)
rpe1_sig_up <- lapply(perturbseq.rpe1, function(x) x$up)

# ------------------------------------------------
# Load perturb-seq data
# ------------------------------------------------

# K562
k562_counts <- read.csv(file.path(PATH, "data/perturb_seq/k562_processed_pb.csv"), row.names = 1)
colnames(k562_counts) <- str_replace_all(colnames(k562_counts), "\\.", "-")

k562_meta <- read.csv(file.path(PATH, "data/perturb_seq/k562_processed_pb_metadata.csv"), row.names = 1)
stopifnot(all.equal(colnames(k562_counts), rownames(k562_meta)))

# RPE1
rpe1_counts <- read.csv(file.path(PATH, "data/perturb_seq/rpe1_processed_pb.csv"), row.names = 1)
colnames(rpe1_counts) <- str_replace_all(colnames(rpe1_counts), "\\.", "-")

rpe1_meta <- read.csv(file.path(PATH, "data/perturb_seq/rpe1_processed_pb_metadata.csv"), row.names = 1)
rownames(rpe1_meta) <- str_replace_all(rownames(rpe1_meta),
                                       "AC118549\\.", "AC118549-")
stopifnot(all.equal(colnames(rpe1_counts), rownames(rpe1_meta)))

# Convert to SummarizedExperiment
k562_se <- SummarizedExperiment(
  assays = list(counts = as.matrix(k562_counts)),
  colData = k562_meta
)

rpe1_se <- SummarizedExperiment(
  assays = list(counts = as.matrix(rpe1_counts)),
  colData = rpe1_meta
)
# ------------------------------------------------
# Load perturbation splits
# ------------------------------------------------

pb_splits <- read.csv(file.path(PATH, "data/sigs/perturb-seq/pb_splits.csv"),
  stringsAsFactors = FALSE
)

# # ------------------------------------------------
# # 1. Control experiment
# # ------------------------------------------------
# 
# k562_ctrl_samples <- rownames(k562_meta)[k562_meta$gene == "non-targeting"]
# rpe1_ctrl_samples <- rownames(rpe1_meta)[rpe1_meta$gene == "non-targeting"]
# 
# k562_ctrl_se <- k562_se[, k562_ctrl_samples]
# rpe1_ctrl_se <- rpe1_se[, rpe1_ctrl_samples]
# 
# # recontextualization
# rpe1_k562_ctrl_gsva  <- projectCor(k562_ctrl_se, rpe1_sig_up, "gsva")
# rpe1_k562_ctrl_eigen <- projectCor(k562_ctrl_se, rpe1_sig_up, "eigen")
# 
# k562_rpe1_ctrl_gsva  <- projectCor(rpe1_ctrl_se, k562_sig_up, "gsva")
# k562_rpe1_ctrl_eigen <- projectCor(rpe1_ctrl_se, k562_sig_up, "eigen")
# 
# saveRDS(rpe1_k562_ctrl_gsva,
#         file.path(PATH,"data/sigs/perturb-seq/projectCor/rpe1_k562_ctrl_gsva.rds"))
# saveRDS(rpe1_k562_ctrl_eigen,
#         file.path(PATH,"data/sigs/perturb-seq/projectCor/rpe1_k562_ctrl_eigen.rds"))
# saveRDS(k562_rpe1_ctrl_gsva,
#         file.path(PATH,"data/sigs/perturb-seq/projectCor/k562_rpe1_ctrl_gsva.rds"))
# saveRDS(k562_rpe1_ctrl_eigen,
#         file.path(PATH,"data/sigs/perturb-seq/projectCor/k562_rpe1_ctrl_eigen.rds"))

# ------------------------------------------------
# 2. 1/10th and 9/10th perturbation experiments
# ------------------------------------------------

foreach(i = c(3,10)) %dopar% {
  
  split_col <- paste0("split_", i)
  message("Processing ", split_col)
  
  # 1/10th genes
  genes_10th <- pb_splits$pb[pb_splits[[split_col]]]
  
  # 9/10th genes
  genes_90th <- pb_splits$pb[!pb_splits[[split_col]]]
  
  # # ------------------------------------------------
  # # 1/10th datasets (control + 1/10th perturbations)
  # # ------------------------------------------------
  # 
  # k562_samples_10th <- rownames(k562_meta)[
  #   k562_meta$gene == "non-targeting" |
  #     k562_meta$gene %in% genes_10th
  # ]
  # 
  # rpe1_samples_10th <- rownames(rpe1_meta)[
  #   rpe1_meta$gene == "non-targeting" |
  #     rpe1_meta$gene %in% genes_10th
  # ]
  # 
  # k562_se_10th <- k562_se[, k562_samples_10th]
  # rpe1_se_10th <- rpe1_se[, rpe1_samples_10th]
  # 
  # rpe1_k562_10th_gsva  <- projectCor(k562_se_10th, rpe1_sig_up, "gsva")
  # rpe1_k562_10th_eigen <- projectCor(k562_se_10th, rpe1_sig_up, "eigen")
  # 
  # k562_rpe1_10th_gsva  <- projectCor(rpe1_se_10th, k562_sig_up, "gsva")
  # k562_rpe1_10th_eigen <- projectCor(rpe1_se_10th, k562_sig_up, "eigen")
  # 
  # saveRDS(rpe1_k562_10th_gsva,
  #         file.path(PATH,paste0("data/sigs/perturb-seq/projectCor/rpe1_k562_10th_gsva_split_",i,".rds")))
  # saveRDS(rpe1_k562_10th_eigen,
  #         file.path(PATH,paste0("data/sigs/perturb-seq/projectCor/rpe1_k562_10th_eigen_split_",i,".rds")))
  # 
  # saveRDS(k562_rpe1_10th_gsva,
  #         file.path(PATH,paste0("data/sigs/perturb-seq/projectCor/k562_rpe1_10th_gsva_split_",i,".rds")))
  # saveRDS(k562_rpe1_10th_eigen,
  #         file.path(PATH,paste0("data/sigs/perturb-seq/projectCor/k562_rpe1_10th_eigen_split_",i,".rds")))
  
  # ------------------------------------------------
  # 9/10th datasets (control + 9/10th perturbations)
  # ------------------------------------------------
  
  k562_samples_90th <- rownames(k562_meta)[
    k562_meta$gene == "non-targeting" |
      k562_meta$gene %in% genes_90th
  ]
  
  rpe1_samples_90th <- rownames(rpe1_meta)[
    rpe1_meta$gene == "non-targeting" |
      rpe1_meta$gene %in% genes_90th
  ]
  
  k562_se_90th <- k562_se[, k562_samples_90th]
  rpe1_se_90th <- rpe1_se[, rpe1_samples_90th]
  
  rpe1_k562_90th_gsva  <- projectCor(k562_se_90th, rpe1_sig_up, "gsva")
  rpe1_k562_90th_eigen <- projectCor(k562_se_90th, rpe1_sig_up, "eigen")
  
  k562_rpe1_90th_gsva  <- projectCor(rpe1_se_90th, k562_sig_up, "gsva")
  k562_rpe1_90th_eigen <- projectCor(rpe1_se_90th, k562_sig_up, "eigen")
  
  saveRDS(rpe1_k562_90th_gsva,
          file.path(PATH,paste0("data/sigs/perturb-seq/projectCor/rpe1_k562_90th_gsva_split_",i,".rds")))
  saveRDS(rpe1_k562_90th_eigen,
          file.path(PATH,paste0("data/sigs/perturb-seq/projectCor/rpe1_k562_90th_eigen_split_",i,".rds")))
  
  saveRDS(k562_rpe1_90th_gsva,
          file.path(PATH,paste0("data/sigs/perturb-seq/projectCor/k562_rpe1_90th_gsva_split_",i,".rds")))
  saveRDS(k562_rpe1_90th_eigen,
          file.path(PATH,paste0("data/sigs/perturb-seq/projectCor/k562_rpe1_90th_eigen_split_",i,".rds")))
}