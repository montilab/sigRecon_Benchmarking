library(tidyverse)
library(doParallel)
library(sigrecon)
library(igraph)
detectCores()
registerDoParallel(15)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
save_path <- file.path(PATH, "data/sigs/drugmatrix/networkprop/")
do_save <- FALSE

# Loading signatures
drugmatrix_sigs <- list(kidney = drugmatrix.kidney,
                        liver = drugmatrix.liver)

# Loading drug splits
drug_splits <- read.csv(file.path(PATH, "data/sigs/drugmatrix/drug_splits.csv"))

# 1. Control Benchmarking
# Loading igraphs
ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/drugmatrix/*control_wgcna.rds"))
contexts <- sub("_control_wgcna\\.rds$", "", basename(ig_paths))
igs <- lapply(ig_paths, readRDS)
names(igs) <- rev(contexts)

recon_dfs <- foreach(context = names(igs), .combine=dplyr::bind_rows) %do% {
  print(context)
  ig <- igs[[context]]
  source_sigs <- lapply(drugmatrix_sigs[[context]], function(x) x$up)
  sig_lengths <- lapply(source_sigs, length) %>% as.numeric
  
  recon_sigs <- network_sig(ig = ig,
                            seeds = source_sigs,
                            sig = "rwr",
                            avg_p = TRUE,
                            bootstrap = TRUE,
                            n_bootstraps = 30,
                            limit = sig_lengths)
  
  saveRDS(recon_sigs, file.path(save_path, paste0(context,"_control.rds")))
}

rm(igs)
print("Done with Control Benchmarking.")

# 1. 1/10th Benchmarking
# Loading igraphs
ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/drugmatrix/*control_1_10th*.rds"))
contexts <- sub("_wgcna.rds$", "", basename(ig_paths))
igs <- lapply(ig_paths, readRDS)
names(igs) <- rev(contexts)

recon_dfs <- foreach(context = names(igs), .combine=dplyr::bind_rows) %dopar% {
  print(context)
  ig <- igs[[context]]
  # Pick appropriate drugs for testing
  split_num <- sub(".*_split_(\\d+)", "\\1", context)
  all_drugs <- drug_splits$drug
  split_col <- paste0("split_", split_num)
  # Get drugs that are FALSE for this split (Network was learned on TRUE)
  drugs_in_split <- drug_splits$drug[drug_splits[[split_col]] == FALSE]
  
  tissue <- sub("_.*", "", context)
  source_sigs <- lapply(drugmatrix_sigs[[tissue]], function(x) x$up)
  shared_drugs <- intersect(drugs_in_split, names(source_sigs))
  source_sigs <- source_sigs[shared_drugs]
  sig_lengths <- lapply(source_sigs, length) %>% as.numeric
  
  recon_sigs <- network_sig(ig = ig,
                            seeds = source_sigs,
                            sig = "rwr",
                            avg_p = TRUE,
                            bootstrap = TRUE,
                            n_bootstraps = 30,
                            limit = sig_lengths)
  
  saveRDS(recon_sigs, file.path(save_path, paste0(context,".rds")))
}

rm(igs)
print("Done with 1/10th Benchmarking.")

# 1. 9/10th Benchmarking
# Loading igraphs
ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/drugmatrix/*control_9_10th*.rds"))
contexts <- sub("_wgcna.rds$", "", basename(ig_paths))
igs <- lapply(ig_paths, readRDS)
names(igs) <- rev(contexts)

recon_dfs <- foreach(context = names(igs), .combine=dplyr::bind_rows) %dopar% {
  print(context)
  ig <- igs[[context]]
  # Pick appropriate drugs for testing
  split_num <- sub(".*_split_(\\d+)", "\\1", context)
  all_drugs <- drug_splits$drug
  split_col <- paste0("split_", split_num)
  # Get drugs that are TRUE for this split (Network was learned on FALSE)
  drugs_in_split <- drug_splits$drug[drug_splits[[split_col]] == TRUE]
  
  tissue <- sub("_.*", "", context)
  source_sigs <- lapply(drugmatrix_sigs[[tissue]], function(x) x$up)
  shared_drugs <- intersect(drugs_in_split, names(source_sigs))
  source_sigs <- source_sigs[shared_drugs]
  sig_lengths <- lapply(source_sigs, length) %>% as.numeric
  
  recon_sigs <- network_sig(ig = ig,
                            seeds = source_sigs,
                            sig = "rwr",
                            avg_p = TRUE,
                            bootstrap = TRUE,
                            n_bootstraps = 30,
                            limit = sig_lengths)
  
  saveRDS(recon_sigs, file.path(save_path, paste0(context,".rds")))
}

rm(igs)
print("Done with 9/10th Benchmarking.")
