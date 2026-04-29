library(tidyverse)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")

tabula_senis_sigs <- readRDS(file.path(PATH, "data/sigs/aging/tabula_muris_senis_sigs_h.rds"))
