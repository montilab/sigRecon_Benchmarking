library(tidyverse)
library(doParallel)
library(sigrecon)
library(igraph)
detectCores()
registerDoParallel(24)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
save_path <- file.path(PATH, "data/sigs/tahoe/networkprop/")
do_save <- FALSE

# Loading signatures
nci_h23_sigs_all <- sigrecon::tahoe.nci_h23

# Loading drug splits
drug_splits <- read.csv(file.path(PATH, "data/sigs/tahoe/drug_splits.csv"))

# # 1. Control Benchmarking
# # Loading igraphs
# ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/tahoe/*control_wgcna.rds"))
# ig_paths <- ig_paths[c(1:5,7:10)]
# cell_lines <- sub("_control_wgcna\\.rds$", "", basename(ig_paths))
# names(ig_paths) <- cell_lines
# 
# 
# # All benchmarking is done on NCI H23 to the other 9 cell lines
# foreach(cell_line = names(ig_paths)) %dopar% {
#   print(cell_line)
#   ig <- readRDS(ig_paths[[cell_line]])
#   source_sigs <- lapply(nci_h23_sigs_all, function(x) x$up)
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
#   saveRDS(recon_sigs, file.path(save_path, paste0(cell_line,"_control.rds")))
#   rm(ig)
#   gc()
# }
# 
# print("Done with Control Benchmarking.")

# ---------------------------
# 1/10th Benchmarking
# ---------------------------

ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/tahoe/*control_1_10th*.rds"))
ig_paths <- ig_paths[-c(51:60)]

cell_lines <- sub("_control_1_10th_split_.*", "", basename(ig_paths))
split_nums <- sub(".*_split_(\\d+)_wgcna\\.rds$", "\\1", basename(ig_paths))

ig_df <- data.frame(
  path = ig_paths,
  cell_line = cell_lines,
  split = split_nums,
  stringsAsFactors = FALSE
)

results <- foreach(cl = unique(ig_df$cell_line)) %:%
  foreach(i = 1:10, .combine = 'c') %dopar% {
    
    split_row <- ig_df %>%
      dplyr::filter(cell_line == cl, split == i)
    
    if(nrow(split_row) == 0) return(NULL)
    
    ig <- readRDS(split_row$path)
    
    split_col <- paste0("split_", i)
    
    drugs_in_split <- drug_splits$drug[drug_splits[[split_col]] == FALSE]
    
    source_sigs <- lapply(nci_h23_sigs_all, function(x) x$up)
    
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
    
    rm(ig)
    gc()
    
    list(list(
      cell_line = cl,
      split = paste0("split_", i),
      sigs = recon_sigs
    ))
  }

# rebuild per-cell-line objects
results_df <- bind_rows(lapply(results, as.data.frame))

for(cl in unique(results_df$cell_line)){
  
  cl_rows <- results[ sapply(results, function(x) x$cell_line == cl) ]
  
  recon_list <- setNames(
    lapply(cl_rows, function(x) x$sigs),
    sapply(cl_rows, function(x) x$split)
  )
  
  saveRDS(
    recon_list,
    file.path(save_path, paste0(cl, "_control_1_10th.rds"))
  )
}

print("Done with 1/10th Benchmarking.")

# ---------------------------
# 9/10th Benchmarking
# ---------------------------

ig_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/tahoe/*control_9_10th*.rds"))
ig_paths <- ig_paths[-c(51:60)]

cell_lines <- sub("_control_9_10th_split_.*", "", basename(ig_paths))
split_nums <- sub(".*_split_(\\d+)_wgcna\\.rds$", "\\1", basename(ig_paths))

ig_df <- data.frame(
  path = ig_paths,
  cell_line = cell_lines,
  split = split_nums,
  stringsAsFactors = FALSE
)

results <- foreach(cl = unique(ig_df$cell_line)) %:%
  foreach(i = 1:10, .combine = 'c') %dopar% {
    
    split_row <- ig_df %>%
      dplyr::filter(cell_line == cl, split == i)
    
    if(nrow(split_row) == 0) return(NULL)
    
    ig <- readRDS(split_row$path)
    
    split_col <- paste0("split_", i)
    
    drugs_in_split <- drug_splits$drug[drug_splits[[split_col]] == TRUE]
    
    source_sigs <- lapply(nci_h23_sigs_all, function(x) x$up)
    
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
    
    rm(ig)
    gc()
    
    list(list(
      cell_line = cl,
      split = paste0("split_", i),
      sigs = recon_sigs
    ))
  }

for(cl in unique(sapply(results, `[[`, "cell_line"))){
  
  cl_rows <- results[sapply(results, function(x) x$cell_line == cl)]
  
  recon_list <- setNames(
    lapply(cl_rows, function(x) x$sigs),
    sapply(cl_rows, function(x) x$split)
  )
  
  saveRDS(
    recon_list,
    file.path(save_path, paste0(cl, "_control_9_10th.rds"))
  )
}

print("Done with 9/10th Benchmarking.")
                        