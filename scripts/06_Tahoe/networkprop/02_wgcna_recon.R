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

# 2. Split benchmarking
run_benchmark <- function(pattern, split_keep, outfile){
  
  ig_paths <- Sys.glob(file.path(PATH, pattern))
  ig_paths <- ig_paths[-c(51:60)]
  
  cell_lines <- sub("_control_.*", "", basename(ig_paths))
  split_nums <- as.integer(sub(".*_split_(\\d+)_wgcna\\.rds$", "\\1", basename(ig_paths)))
  
  ig_df <- data.frame(
    path = ig_paths,
    cell_line = cell_lines,
    split = split_nums,
    stringsAsFactors = FALSE
  )
  
  source_sigs_full <- lapply(nci_h23_sigs_all, function(x) x$up)
  
  results <- foreach(
    j = seq_len(nrow(ig_df)),
    .packages=c("sigrecon","igraph"),
    .combine=function(a,b){
      for(cl in names(b)){
        if(is.null(a[[cl]])) a[[cl]] <- list()
        a[[cl]] <- c(a[[cl]], b[[cl]])
      }
      a
    },
    .init=list()
  ) %dopar% {
    
    cl <- ig_df$cell_line[j]
    split_i <- ig_df$split[j]
    path <- ig_df$path[j]
    
    ig <- readRDS(path)
    
    split_col <- paste0("split_", split_i)
    
    drugs_in_split <- drug_splits$drug[drug_splits[[split_col]] == split_keep]
    
    shared_drugs <- intersect(drugs_in_split, names(source_sigs_full))
    source_sigs <- source_sigs_full[shared_drugs]
    
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
    
    list(
      setNames(
        list(setNames(list(recon_sigs), paste0("split_", split_i))),
        cl
      )
    )
  }
  
  saveRDS(results, file.path(save_path, outfile))
}

# run experiments
run_benchmark(
  "data/wgcna_networks/tahoe/*control_1_10th*.rds",
  FALSE,
  "control_1_10th.rds"
)

run_benchmark(
  "data/wgcna_networks/tahoe/*control_9_10th*.rds",
  TRUE,
  "control_9_10th.rds"
)
