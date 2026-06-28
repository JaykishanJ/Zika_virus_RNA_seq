# ==============================================================================
# Meta-Analysis Script (Meta_Analysis.R)
# Combines GSE146423 and GSE265922 raw counts into a unified DESeq2 model
# ==============================================================================

library(here)
library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)

# Ensure we're in the right place relative to the project root
# If ran interactively, ensure working directory is set to Zika_virus_wetlab

cat("============================================================\n")
cat("ZIKV Meta-Analysis (GSE146423 + GSE265922)\n")
cat("============================================================\n\n")

# Create output directories
out_dir <- here::here("Meta_Analysis", "Results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "QC"), showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), showWarnings = FALSE)

# 1. Load Raw Count Matrices
cat("[1] Loading count matrices...\n")
counts_146423_file <- here::here("GSE146423", "Results", "Counts", "GSE146423_RAW_NCBI_Count_Matrix.csv")
counts_265922_file <- here::here("GSE265922", "Results", "Counts", "GSE265922_RAW_STAR_Unstranded_Count_Matrix.csv")

if (!file.exists(counts_146423_file)) {
  stop("GSE146423 count matrix not found. Did you run DEG_GSE146423.R?")
}
if (!file.exists(counts_265922_file)) {
  stop("GSE265922 count matrix not found. Did you run DEG_GES265922.R?")
}

counts_146423 <- read.csv(counts_146423_file, row.names = 1, check.names = FALSE)
counts_265922 <- read.csv(counts_265922_file, row.names = 1, check.names = FALSE)

# Convert to data frames with a gene_id column for merging
df_146423 <- rownames_to_column(counts_146423, "gene_id")
df_265922 <- rownames_to_column(counts_265922, "gene_id")

# 2. Find common genes and merge
cat("[2] Merging count matrices on common ENSEMBL/Gene IDs...\n")
# Note: GSE146423 might use different gene identifiers depending on the GEO file.
# We will inner join on whatever identifier is in the rownames. 
# Usually, GSE265922 is ENSEMBL (e.g. ENSG...). GSE146423 from NCBI might be GeneID/Symbol.
# If they don't overlap well, a mapping step would be required. 
# We perform inner join assuming both are ENSEMBL or GeneID.

combined_counts <- inner_join(df_146423, df_265922, by = "gene_id")

if (nrow(combined_counts) == 0) {
  warning("No overlapping gene IDs found. GSE146423 and GSE265922 might use different ID formats (e.g. Symbol vs ENSEMBL). Gene ID mapping will be required for a full meta-analysis.")
} else {
  cat("Found", nrow(combined_counts), "overlapping genes.\n")
  
  count_mat <- as.matrix(column_to_rownames(combined_counts, "gene_id"))
  
  # 3. Create combined Phenotype table
  cat("[3] Creating combined phenotype...\n")
  # GSE146423: 3 Control, 3 ZIKV
  # GSE265922: 3 Mock, 3 ZIKV
  
  samples_146423 <- colnames(counts_146423)
  samples_265922 <- colnames(counts_265922)
  
  pheno_146423 <- data.frame(
    Sample = samples_146423,
    Dataset = "GSE146423",
    Condition = ifelse(grepl("ZIKV", samples_146423), "ZIKV", "Control") # Adjust based on sample naming
  )
  
  pheno_265922 <- data.frame(
    Sample = samples_265922,
    Dataset = "GSE265922",
    Condition = ifelse(grepl("^Z", samples_265922), "ZIKV", "Control") # Assuming Z1, Z2, Z3
  )
  
  pheno_combined <- bind_rows(pheno_146423, pheno_265922)
  rownames(pheno_combined) <- pheno_combined$Sample
  
  # Ensure factors
  pheno_combined$Dataset <- factor(pheno_combined$Dataset)
  pheno_combined$Condition <- factor(pheno_combined$Condition, levels = c("Control", "ZIKV"))
  
  # Ensure column order matches
  count_mat <- count_mat[, rownames(pheno_combined)]
  
  # 4. Run DESeq2 with Batch Correction
  cat("[4] Running DESeq2 with ~ Dataset + Condition...\n")
  dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = pheno_combined, design = ~ Dataset + Condition)
  
  # Filter low count genes
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("Condition", "ZIKV", "Control"), alpha = 0.05)
  
  res_shrunk <- lfcShrink(dds, coef = "Condition_ZIKV_vs_Control", type = "apeglm")
  
  # Save Results
  write.csv(as.data.frame(res_shrunk), file.path(out_dir, "Tables", "Meta_Analysis_DESeq2_ZIKV_vs_Control.csv"))
  
  cat("Meta-Analysis successfully generated unified results!\n")
}
