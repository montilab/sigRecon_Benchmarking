library(tidyverse)
library(limma)
library(Biobase)
options(box.path=file.path(Sys.getenv("CBMGIT"), "MLscripts"))

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/drugmatrix/")
do_save <- TRUE

# 0. Load Expression Sets
liver_eset <- readRDS(file=file.path(DATA_PATH, "liver.rds"))
kidney_eset <- readRDS(file=file.path(DATA_PATH, "kidney.rds"))
print(dim(liver_eset))
print(dim(kidney_eset))

# 1. Remove lowly expressed genes
box::use(R/rm_low_rnaseq_counts)
liver_eset_filtered <- rm_low_rnaseq_counts$rm_low_rnaseq_counts(liver_eset, min_samples = 3)
kidney_eset_filtered <- rm_low_rnaseq_counts$rm_low_rnaseq_counts(kidney_eset, min_samples = 3)
print(dim(liver_eset_filtered))
print(dim(kidney_eset_filtered))

# No genes are filtered

# Plot Variables
plot_variables <- c("rna extraction date:ch1", "time:ch1")

# ============================================================================
# 2. Select Top HVGs and Perform PCA
# ============================================================================

#' Select top HVGs by variance and perform PCA
#' @param eset ExpressionSet object
#' @param n_hvgs Number of highly variable genes to select
#' @return List containing PCA results and variance explained
select_hvgs_and_pca <- function(eset, n_hvgs = 3000) {
  
  exprs_mat <- exprs(eset)
  
  # Calculate variance for each gene
  gene_vars <- apply(exprs_mat, 1, var)
  
  # Select top N HVGs
  n_hvgs <- min(n_hvgs, nrow(exprs_mat))  # Don't exceed total genes
  top_hvg_indices <- order(gene_vars, decreasing = TRUE)[1:n_hvgs]
  hvg_names <- rownames(exprs_mat)[top_hvg_indices]
  
  # Subset to HVGs
  exprs_hvg <- exprs_mat[top_hvg_indices, ]
  
  cat("Selected", n_hvgs, "HVGs out of", nrow(exprs_mat), "genes\n")
  cat("Variance range of HVGs:", 
      round(min(gene_vars[top_hvg_indices]), 2), "to",
      round(max(gene_vars[top_hvg_indices]), 2), "\n")
  
  # Transpose for PCA (samples as rows)
  exprs_t <- t(exprs_hvg)
  
  # Perform PCA
  pca_result <- prcomp(exprs_t, center = TRUE, scale. = TRUE)
  
  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, ] * 100
  
  cat("Variance explained by PC1-PC5:",
      paste(round(var_explained[1:5], 1), "%", collapse = ", "), "\n\n")
  
  return(list(
    pca = pca_result,
    var_explained = var_explained,
    hvg_names = hvg_names,
    hvg_indices = top_hvg_indices
  ))
}

#' Create PCA plots colored by batch variables
#' @param eset ExpressionSet object
#' @param pca_result PCA results from select_hvgs_and_pca
#' @param batch_vars Character vector of batch variable names
#' @param dataset_name Name for plot titles
#' @return List of ggplot objects
create_pca_plots <- function(eset, pca_result, batch_vars, dataset_name) {
  
  # Extract PC coordinates
  pc_coords <- as.data.frame(pca_result$pca$x)
  
  # Add metadata
  pdata <- pData(eset)
  pc_data <- cbind(pc_coords, pdata)
  
  # Variance explained for axis labels
  var_explained <- pca_result$var_explained
  
  plot_list <- list()
  
  # Create plot for each batch variable
  for (batch_var in batch_vars) {
    
    # Check if variable exists
    if (!batch_var %in% colnames(pdata)) {
      warning("Batch variable '", batch_var, "' not found in metadata. Skipping.")
      next
    }
    
    # Check for NA values
    if (all(is.na(pc_data[[batch_var]]))) {
      warning("All values are NA for '", batch_var, "'. Skipping.")
      next
    }
    
    # Convert to factor if character or factor
    if (is.character(pc_data[[batch_var]]) || is.factor(pc_data[[batch_var]])) {
      pc_data[[batch_var]] <- as.factor(pc_data[[batch_var]])
      
      # If too many levels, show warning
      n_levels <- length(unique(pc_data[[batch_var]]))
      if (n_levels > 20) {
        warning("Batch variable '", batch_var, "' has ", n_levels, 
                " levels. Plot may be cluttered.")
      }
    }
    
    # Create plot
    p <- ggplot(pc_data, aes(x = PC1, y = PC2, color = .data[[batch_var]])) +
      geom_point(size = 2.5, alpha = 0.7) +
      labs(
        title = dataset_name,
        subtitle = paste("Colored by:", batch_var),
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
        color = gsub(":ch1", "", batch_var)  # Clean up label
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5)
      )
    
    # Adjust legend based on variable type
    if (is.numeric(pc_data[[batch_var]])) {
      p <- p + scale_color_viridis_c(option = "plasma")
    } else {
      # Use colorblind-friendly palette for categorical
      n_colors <- length(unique(pc_data[[batch_var]][!is.na(pc_data[[batch_var]])]))
      if (n_colors <= 8) {
        p <- p + scale_color_brewer(palette = "Set2", na.value = "grey80")
      } else if (n_colors <= 12) {
        p <- p + scale_color_brewer(palette = "Paired", na.value = "grey80")
      } else {
        p <- p + scale_color_viridis_d(option = "turbo", na.value = "grey80")
      }
    }
    
    plot_list[[batch_var]] <- p
  }
  
  return(plot_list)
}

# ============================================================================
# 3. Run PCA Analysis
# ============================================================================
# Liver PCA
liver_pca <- select_hvgs_and_pca(liver_eset_filtered, n_hvgs = 3000)
liver_plots <- create_pca_plots(
  liver_eset_filtered, 
  liver_pca, 
  plot_variables, 
  "Liver (Top 3000 HVGs)"
)
print(liver_plots)

# Kidney PCA
kidney_pca <- select_hvgs_and_pca(kidney_eset_filtered, n_hvgs = 3000)
kidney_plots <- create_pca_plots(
  kidney_eset_filtered, 
  kidney_pca, 
  plot_variables, 
  "Kidney (Top 3000 HVGs)"
)
print(kidney_plots)
