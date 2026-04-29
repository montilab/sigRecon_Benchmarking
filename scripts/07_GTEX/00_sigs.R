library(tidyverse)
el_aging_sig_df <- read.csv("/restricted/projectnb/brcameta/projects/sig_recon/data/sigs/aging/Combined_Aging_LLFS_ILO_Mortality_LLFS.csv", row.names = 1)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")

ilo_aging_sig <- el_aging_sig_df %>% 
  dplyr::filter(age_sig_ilo) %>%
  dplyr::filter(age_dir_ilo == "Inc.") %>%
  dplyr::arrange(age_qval_ilo) %>%
  dplyr::pull(gene_id)
ilo_aging_sig_list <- list(up = ilo_aging_sig[1:100], up_full = ilo_aging_sig)

llfs_aging_sig <- el_aging_sig_df %>% 
  dplyr::filter(age_sig_llfs) %>%
  dplyr::filter(age_dir_llfs == "Inc.") %>%
  dplyr::arrange(age_qval_llfs) %>%
  dplyr::pull(gene_id)

llfs_aging_sig_list <- list(up = llfs_aging_sig[1:100], up_full = llfs_aging_sig)

saveRDS(ilo_aging_sig_list, file.path(PATH, "data/sigs/aging/ilo_aging_sig.rds"))
saveRDS(llfs_aging_sig_list, file.path(PATH, "data/sigs/aging/llfs_aging_sig.rds"))
