library(tidyseurat)
library(Seurat)
library(sceasy)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(reticulate)
reticulate::use_condaenv("r-sceasy")

DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/srivatsan_2019")
PROJ_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/stack")

drug_splits <- read.csv("/restricted/projectnb/brcameta/projects/sig_recon/data/sigs/sciplex/drug_splits.csv")

a549 <- readRDS(file.path(DATA_PATH, "a549.rds"))
k562 <- readRDS(file.path(DATA_PATH, "k562.rds"))
mcf7 <- readRDS(file.path(DATA_PATH, "mcf7.rds"))

a549_ctrl <- a549 %>% dplyr::filter(product_name == "Vehicle")
a549_base <- a549 %>% dplyr::filter(dose == 10000 | product_name == "Vehicle")
k562_ctrl <- k562 %>% dplyr::filter(product_name == "Vehicle")
k562_base <- k562 %>% dplyr::filter(dose == 10000 | product_name == "Vehicle")
mcf7_ctrl <- mcf7 %>% dplyr::filter(product_name == "Vehicle")
mcf7_base <- mcf7 %>% dplyr::filter(dose == 10000 | product_name == "Vehicle")

sceasy::convertFormat(a549_ctrl, from="seurat", to="anndata",
                      main_layer = "data", # main_layer gets copied to adata.X, in this case data are counts
                      outFile=file.path(PROJ_PATH, "sciplex/a549_ctrl.h5ad"))
sceasy::convertFormat(a549_base, from="seurat", to="anndata",
                      main_layer = "data", # main_layer gets copied to adata.X, in this case data are counts
                      outFile=file.path(PROJ_PATH, "sciplex/a549_all.h5ad"))
sceasy::convertFormat(k562_ctrl, from="seurat", to="anndata",
                      main_layer = "data", # main_layer gets copied to adata.X, in this case data are counts
                      outFile=file.path(PROJ_PATH, "sciplex/k562_ctrl.h5ad"))
sceasy::convertFormat(k562_base, from="seurat", to="anndata",
                      main_layer = "data", # main_layer gets copied to adata.X, in this case data are counts
                      outFile=file.path(PROJ_PATH, "sciplex/k562_all.h5ad"))
sceasy::convertFormat(mcf7_ctrl, from="seurat", to="anndata",
                      main_layer = "data", # main_layer gets copied to adata.X, in this case data are counts
                      outFile=file.path(PROJ_PATH, "sciplex/mcf7_ctrl.h5ad"))
sceasy::convertFormat(mcf7_base, from="seurat", to="anndata",
                      main_layer = "data", # main_layer gets copied to adata.X, in this case data are counts
                      outFile=file.path(PROJ_PATH, "sciplex/mcf7_all.h5ad"))