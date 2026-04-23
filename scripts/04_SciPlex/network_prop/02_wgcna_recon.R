library(tidyverse)
library(doParallel)
library(sigrecon)
library(igraph)
detectCores()
registerDoParallel(15)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
save_path <- file.path(PATH, "data/sigs/sciplex/networkprop")
do_save <- FALSE

# Loading signatures
sciplex_sigs <- list(a549 = sciplex.a549,
                     k562 = sciplex.k562,
                     mcf7 = sciplex.mcf7)

# Loading drug splits
drug_splits <- read.csv(file.path(PATH, "data/sigs/sciplex/drug_splits.csv"))

# 1. Control Benchmarking
# Loading igraphs
ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/sciplex/*control_wgcna.rds"))
contexts <- sub("_control_wgcna\\.rds$", "", basename(ig_paths))
# igs <- lapply(ig_paths, readRDS)
# names(igs) <- contexts

comparisons <- expand.grid(contexts, contexts) %>%
  dplyr::rename(source = Var1, target = Var2) %>%
  dplyr::filter(source != target)
# 
# foreach(i = 1:nrow(comparisons)) %dopar% {
#   source <- comparisons[i, "source"]
#   target <- comparisons[i, "target"]
#   
#   ig <- igs[[target]]
#   source_sigs <- lapply(sciplex_sigs[[source]], function(x) x$up)
#   sig_lengths <- lapply(source_sigs, length) %>% as.numeric
# 
#   recon_sigs <- network_sig(ig = ig,
#                             seeds = source_sigs,
#                             sig = "rwr",
#                             avg_p = TRUE,
#                             bootstrap = TRUE,
#                             n_bootstraps = 30,
#                             limit = sig_lengths)
#   saveRDS(recon_sigs, file.path(save_path, paste0(source,"_",target,"_control.rds")))
# }
# rm(igs)
# print("Done with Control Benchmarking.")

# 10th benchmarking
ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/sciplex/*control_1_10th*.rds"))
split_names <- sub("_wgcna.rds$", "", basename(ig_paths))
names(ig_paths) <- split_names

foreach(i = 1:nrow(comparisons)) %dopar% {
  
  source <- comparisons[i, "source"]
  target <- comparisons[i, "target"]
  
  message("Processing ", source, " -> ", target)
  
  # graphs for this target cell line
  target_graphs <- ig_paths[grep(paste0("^", target, "_"), names(ig_paths))]
  
  recon_list <- list()
  
  for(split_name in names(target_graphs)) {
    
    ig <- readRDS(target_graphs[[split_name]])
    
    split_num <- sub(".*_split_(\\d+)", "\\1", split_name)
    split_col <- paste0("split_", split_num)
    
    drugs_in_split <- drug_splits$drug[drug_splits[[split_col]] == FALSE]
    
    source_sigs <- lapply(sciplex_sigs[[source]], function(x) x$up)
    
    shared_drugs <- intersect(drugs_in_split, names(source_sigs))
    source_sigs <- source_sigs[shared_drugs]
    
    sig_lengths <- lengths(source_sigs)
    
    recon_sigs <- network_sig(
      ig = ig,
      seeds = source_sigs,
      sig = "rwr",
      avg_p = TRUE,
      bootstrap = TRUE,
      n_bootstraps = 30,
      limit = sig_lengths
    )
    
    recon_list[[paste0("split_", split_num)]] <- recon_sigs
    
    rm(ig)
    gc()
  }
  
  saveRDS(
    recon_list,
    file.path(save_path, paste0(source, "_", target, "_1_10th.rds"))
  )
}

print("Done with 1/10th Benchmarking.")

ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/sciplex/*control_9_10th*.rds"))
split_names <- sub("_wgcna.rds$", "", basename(ig_paths))
names(ig_paths) <- split_names

foreach(i = 1:nrow(comparisons)) %dopar% {
  
  source <- comparisons[i, "source"]
  target <- comparisons[i, "target"]
  
  message("Processing ", source, " -> ", target)
  
  target_graphs <- ig_paths[grep(paste0("^", target, "_"), names(ig_paths))]
  
  recon_list <- list()
  
  for(split_name in names(target_graphs)) {
    
    ig <- readRDS(target_graphs[[split_name]])
    
    split_num <- sub(".*_split_(\\d+)", "\\1", split_name)
    split_col <- paste0("split_", split_num)
    
    drugs_in_split <- drug_splits$drug[drug_splits[[split_col]] == TRUE]
    
    source_sigs <- lapply(sciplex_sigs[[source]], function(x) x$up)
    
    shared_drugs <- intersect(drugs_in_split, names(source_sigs))
    source_sigs <- source_sigs[shared_drugs]
    
    sig_lengths <- lengths(source_sigs)
    
    recon_sigs <- network_sig(
      ig = ig,
      seeds = source_sigs,
      sig = "rwr",
      avg_p = TRUE,
      bootstrap = TRUE,
      n_bootstraps = 30,
      limit = sig_lengths
    )
    
    recon_list[[paste0("split_", split_num)]] <- recon_sigs
    
    rm(ig)
    gc()
  }
  
  saveRDS(
    recon_list,
    file.path(save_path, paste0(source, "_", target, "_9_10th.rds"))
  )
}

print("Done with 9/10th Benchmarking.")
