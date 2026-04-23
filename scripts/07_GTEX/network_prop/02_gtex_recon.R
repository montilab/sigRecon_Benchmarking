library(tidyverse)
library(doParallel)
library(sigrecon)
# detectCores()
registerDoParallel(35)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
do_save <- TRUE

# Loading Signatures and Networks
load(file.path(PATH, "sigrecon/data/gtex.aging.rda"))
blood_sig <- gtex.aging["Whole_Blood"]
names(blood_sig) <- "Aging"
blood_ig <- readRDS(file.path(PATH, "data/wgcna_networks/gtex/Whole_Blood.rds"))

all_net_paths <- Sys.glob(file.path(PATH, "data/wgcna_networks/gtex/*.rds"))
non_blood_net_paths <- all_net_paths[str_detect(all_net_paths, "Whole_Blood", negate=TRUE)]

# 1. Baseline: Without Recontextualization
no_recon <- foreach(filepath=non_blood_net_paths, .combine=dplyr::bind_rows) %dopar% {
  tissue_name <- basename(filepath) |> str_remove(pattern = "\\.rds$")
  dest_ig <- readRDS(filepath)
  dest_sig <- gtex.aging[tissue_name]
  names(dest_sig) <- "Aging"
  sig_length <- length(blood_sig[["up"]])

  recon_eval_df(ig = dest_ig,
                seed_name = paste0("Whole_Blood", "_", tissue_name),
                source_sigs = blood_sig,
                dest_sigs = dest_sig,
                restart=-1,
                recon=FALSE,
                use_weights=TRUE,
                weights.pwr=1,
                limit = sig_length)

}
print("Done with Baseline")
saveRDS(no_recon, file=file.path(PATH, "data/eval/gtex/no_recon_wgcna_power.rds"))

# 2. With Recontextualization
restart_vals <- c(1,seq(1e-1,1e-4, length.out = 5))

recon_dfs <- foreach(rw_p = restart_vals, .combine=dplyr::bind_rows) %:%
  foreach(filepath=non_blood_net_paths) %dopar% {
    tissue_name <- basename(filepath) |> str_remove(pattern = "\\.rds$")
    print(paste(rw_p, tissue_name))
    dest_ig <- readRDS(filepath)
    dest_sig <- gtex.aging[tissue_name]
    names(dest_sig) <- "Aging"
    sig_length <- length(blood_sig[["up"]])

    recon_eval_df(ig = dest_ig,
                  seed_name = paste0("Whole_Blood", "_", tissue_name),
                  source_sigs = blood_sig,
                  dest_sigs = dest_sig,
                  restart=rw_p,
                  avg_p = FALSE,
                  bootstrap = TRUE,
                  n_bootstraps = 250,
                  recon=TRUE,
                  use_weights=TRUE,
                  weights.pwr = 1,
                  limit = sig_length)
}
print("Done with single restart values.")
saveRDS(recon_dfs, file=file.path(PATH, "data/eval/gtex/recon_wgcna_power_bootstrap.rds"))

# 3. Averaging restart_values
restart_vals <- c(1,seq(1e-1,1e-4, length.out = 5))

recon_dfs <- foreach(filepath=non_blood_net_paths, .combine=dplyr::bind_rows) %dopar% {
    tissue_name <- basename(filepath) |> str_remove(pattern = "\\.rds$")
    dest_ig <- readRDS(filepath)
    dest_sig <- gtex.aging[tissue_name]
    names(dest_sig) <- "Aging"
    sig_length <- length(blood_sig[["up"]])

    recon_eval_df(ig = dest_ig,
                  seed_name = paste0("Whole_Blood", "_", tissue_name),
                  source_sigs = blood_sig,
                  dest_sigs = dest_sig,
                  avg_p = TRUE,
                  bootstrap = TRUE,
                  n_bootstraps = 250,
                  recon=TRUE,
                  save=TRUE,
                  save_path = file.path(PATH, "data/sigs/gtex/recon/wgcna"),
                  use_weights=TRUE,
                  weights.pwr = 1,
                  limit = sig_length)
  }
print("Done with averaging restart values.")
saveRDS(recon_dfs, file=file.path(PATH, "data/eval/gtex/recon_avg_p_wgcna_power_bootstrap.rds"))
