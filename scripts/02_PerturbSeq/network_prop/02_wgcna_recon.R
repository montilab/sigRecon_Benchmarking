library(tidyverse)
library(doParallel)
library(sigrecon)
library(igraph)
detectCores()
registerDoParallel(20)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
save_path <- file.path(PATH, "data/sigs/perturb-seq/networkprop/")
do_save <- FALSE

# Loading signatures
perturbseq_sigs <- list(rpe1 = perturbseq.rpe1,
                     k562 = perturbseq.k562)

# Loading pb splits
pb_splits <- read.csv(file.path(PATH, "data/sigs/perturb-seq/pb_splits.csv"))

# # # 1. Control Benchmarking
# # Loading igraphs
# ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/perturb-seq/*control_wgcna.rds"))
# contexts <- sub("_control_wgcna\\.rds$", "", basename(ig_paths))
# igs <- lapply(ig_paths, readRDS)
# names(igs) <- rev(contexts)
# 
# recon_dfs <- foreach(context = names(igs), .combine=dplyr::bind_rows) %do% {
#   print(context)
#   ig <- igs[[context]]
#   source_sigs <- lapply(perturbseq_sigs[[context]], function(x) x$up)
#   sig_lengths <- lapply(source_sigs, length) %>% as.numeric
# 
#   recon_sigs <- network_sig(ig = ig,
#                             seeds = source_sigs,
#                             sig = "rwr",
#                             avg_p = TRUE,
#                             bootstrap = TRUE,
#                             n_bootstraps = 30,
#                             limit = sig_lengths)
# 
#   saveRDS(recon_sigs, file.path(save_path, paste0(context,"_control.rds")))
# }
# 
# rm(igs)
# print("Done with Control Benchmarking.")

# 1. 1/10th Benchmarking
# Loading igraphs
ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/perturb-seq/*control_1_10th*.rds"))
contexts <- sub("_wgcna.rds$", "", basename(ig_paths))
names(ig_paths) <- rev(contexts)

recon_dfs <- foreach(context = names(ig_paths), .combine=dplyr::bind_rows) %dopar% {
  print(context)
  ig <- readRDS(ig_paths[[context]])
  # Pick appropriate pbs for testing
  split_num <- sub(".*_split_(\\d+)", "\\1", context)
  all_pbs <- pb_splits$pb
  split_col <- paste0("split_", split_num)
  # Get pbs that are FALSE for this split (Network was learned on TRUE)
  pbs_in_split <- pb_splits$pb[pb_splits[[split_col]] == FALSE]
  
  tissue <- sub("_.*", "", context)
  source_sigs <- lapply(perturbseq_sigs[[tissue]], function(x) x$up)
  shared_pbs <- intersect(pbs_in_split, names(source_sigs))
  source_sigs <- source_sigs[shared_pbs]
  sig_lengths <- lapply(source_sigs, length) %>% as.numeric
  
  recon_sigs <- network_sig(ig = ig,
                            seeds = source_sigs,
                            sig = "rwr",
                            avg_p = TRUE,
                            bootstrap = TRUE,
                            n_bootstraps = 30,
                            limit = sig_lengths)
  
  saveRDS(recon_sigs, file.path(save_path, paste0(context,"._10th.rds")))
  rm(ig)
  gc()
}
print("Done with 1/10th Benchmarking.")

# 1. 9/10th Benchmarking
# Loading igraphs
ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/perturb-seq/*control_9_10th*.rds"))
contexts <- sub("_wgcna.rds$", "", basename(ig_paths))
names(ig_paths) <- rev(contexts)

recon_dfs <- foreach(context = names(ig_paths), .combine=dplyr::bind_rows) %dopar% {
  print(context)
  ig <- readRDS(ig_paths[[context]])
  # Pick appropriate pbs for testing
  split_num <- sub(".*_split_(\\d+)", "\\1", context)
  all_pbs <- pb_splits$pb
  split_col <- paste0("split_", split_num)
  # Get pbs that are TRUE for this split (Network was learned on FALSE)
  pbs_in_split <- pb_splits$pb[pb_splits[[split_col]] == TRUE]
  
  tissue <- sub("_.*", "", context)
  source_sigs <- lapply(perturbseq_sigs[[tissue]], function(x) x$up)
  shared_pbs <- intersect(pbs_in_split, names(source_sigs))
  source_sigs <- source_sigs[shared_pbs]
  sig_lengths <- lapply(source_sigs, length) %>% as.numeric
  
  recon_sigs <- network_sig(ig = ig,
                            seeds = source_sigs,
                            sig = "rwr",
                            avg_p = TRUE,
                            bootstrap = TRUE,
                            n_bootstraps = 30,
                            limit = sig_lengths)
  
  saveRDS(recon_sigs, file.path(save_path, paste0(context,"_90th.rds")))
  rm(ig)
  gc()
}
print("Done with 9/10th Benchmarking.")
