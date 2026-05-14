library(tidyverse)
library(readxl)
library(babelgene)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
el_aging_sig_df <- read.csv(file.path(PATH, "data/sigs/aging/Combined_Aging_LLFS_ILO_Mortality_LLFS.csv", row.names = 1))

## Source signatures
# tabulasenis_sigs <- readRDS("/restricted/projectnb/agedisease/projects/challenge2025/data/Tabula_senis/tms_combined_dir/TabulaSenis_omic_collection_by_tissue.rds")
regeneron_aging_sig <- read_excel(file.path(PATH,"data/sigs/aging/ACEL-25-e70394-s014.xlsx"),sheet = "PBMC")
regeneron_aging_sig$score <- regeneron_aging_sig$Log2FC*(-log10(regeneron_aging_sig$padj))
up_sig <- regeneron_aging_sig %>% 
  dplyr::arrange(desc(score)) %>%
  dplyr::slice(1:100) %>% 
  dplyr::pull(gene)
up_full_sig <- regeneron_aging_sig %>% 
  dplyr::arrange(desc(score)) %>%
  dplyr::pull(gene)
regeneron_sig <- list(up = up_sig, up_full = up_full_sig)
regeneron_sig_h <- lapply(regeneron_sig, function(x) babelgene::orthologs(genes = x, species = "mouse", human = FALSE)$human_ensembl)
saveRDS(regeneron_sig, file.path(PATH, "data/sigs/aging/regeneron_sig_m.rds"))
saveRDS(regeneron_sig_h, file.path(PATH, "data/sigs/aging/regeneron_sig_h.rds"))
## Ground truth signature
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
