library(tidyverse)
library(sigrecon)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")

# Loading Sigs
regeneron_sig_h <- readRDS(file.path(PATH, "data/sigs/aging/regeneron_sig_h.rds"))
ilo_aging_sig <- readRDS(file.path(PATH, "data/sigs/aging/ilo_aging_sig.rds"))
llfs_aging_sig <- readRDS(file.path(PATH, "data/sigs/aging/llfs_aging_sig.rds"))
gtex_aging_sig <- list(aging = gsub("\\.\\d+$", "", gtex.blood$up))
regeneron_up <- list(aging = regeneron_sig_h$up)
ilo_aging_sig <- list(aging = ilo_aging_sig)
llfs_aging_sig <- list(aging = llfs_aging_sig)

# No-change
sig_eval_table(source_sigs = regeneron_up,pred_sigs = regeneron_up,true_sigs = llfs_aging_sig)
sig_eval_table(source_sigs = regeneron_up,pred_sigs = regeneron_up,true_sigs = ilo_aging_sig)

# Mean (just using GTEX)
sig_eval_table(source_sigs = regeneron_up,pred_sigs = gtex_aging_sig,true_sigs = llfs_aging_sig)
sig_eval_table(source_sigs = regeneron_up,pred_sigs = gtex_aging_sig,true_sigs = ilo_aging_sig)

# vennr::vennr(list(LLFS = llfs_aging_sig$aging$up, ILO = ilo_aging_sig$aging$up))
