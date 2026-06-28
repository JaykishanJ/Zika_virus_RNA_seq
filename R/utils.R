# ==============================================================================
# Helper Utilities (utils.R)
# Contains general data manipulation and installation functions
# ==============================================================================

# Function to install missing packages
install_if_missing <- function(pkgs, bioc = FALSE) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing missing package:", pkg, "\n")
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        }
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
  }
}

# Function to clean GSEA data frames
clean_gsea_df <- function(gsea_df) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0) {
    return(gsea_df)
  }
  
  required_cols <- c("ID", "pvalue", "p.adjust", "NES")
  missing_cols <- setdiff(required_cols, colnames(gsea_df))
  
  if (length(missing_cols) > 0) {
    return(gsea_df)
  }
  
  gsea_df %>%
    dplyr::filter(
      !is.na(ID),
      ID != "NA",
      is.finite(pvalue),
      is.finite(p.adjust),
      is.finite(NES)
    )
}
