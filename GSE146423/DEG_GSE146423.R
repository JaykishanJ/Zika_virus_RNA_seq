library(here)
setwd(here::here("GSE146423"))

# ==============================================================================
# 0. CLEAN SESSION AND USER SETTINGS
# ==============================================================================

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
set.seed(123)


GSE_ID <- "GSE146423"

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0

MIN_COUNT   <- 10
MIN_SAMPLES <- 3

TOP_N_MAIN_HEATMAP <- 30
TOP_N_SUPP_HEATMAP <- 50

dir.create("Results", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Raw", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Counts", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/QC", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Plots", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Tables", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Annotation", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Enrichment", showWarnings = FALSE, recursive = TRUE)

cat("\n============================================================\n")
cat("Raw RNA-seq DEG Analysis:", GSE_ID, "\n")
cat("Working directory:", getwd(), "\n")
cat("Comparison: ZIKV infected vs Control\n")
cat("Count type: NCBI/GEO-generated raw gene counts\n")
cat("============================================================\n\n")


# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================

cat("[1] Loading required libraries...\n")

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "scales",
  "RColorBrewer",
  "viridis",
  "forcats",
  "BiocManager"
)

bioc_pkgs <- c(
  "DESeq2",
  "apeglm",
  "EnhancedVolcano",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot"
)

source(here::here("R", "utils.R"))

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(scales)
  library(RColorBrewer)
  library(viridis)
  library(forcats)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
})

cat("✓ Libraries loaded successfully.\n\n")


# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

cat("[2] Defining helper functions...\n")

source(here::here("R", "plot_functions.R"))

cat("✓ Helper functions ready.\n\n")


# ==============================================================================
# 3. DEFINE CORRECT PHENOTYPE / METADATA
# ==============================================================================

cat("[3] Building correct sample phenotype table...\n")

phenotype <- data.frame(
  geo_accession = c(
    "GSM4385296",
    "GSM4385297",
    "GSM4385298",
    "GSM4385299",
    "GSM4385300",
    "GSM4385301"
  ),
  SampleID = c(
    "C1",
    "C2",
    "C3",
    "ZIKV1",
    "ZIKV2",
    "ZIKV3"
  ),
  Condition = c(
    "Control",
    "Control",
    "Control",
    "ZIKV",
    "ZIKV",
    "ZIKV"
  ),
  Replicate = c(
    "rep1",
    "rep2",
    "rep3",
    "rep1",
    "rep2",
    "rep3"
  ),
  Cell_Line = "A549",
  Treatment = c(
    "A549 Control",
    "A549 Control",
    "A549 Control",
    "A549 infected with Zika virus",
    "A549 infected with Zika virus",
    "A549 infected with Zika virus"
  ),
  stringsAsFactors = FALSE
)

phenotype$Condition <- factor(
  phenotype$Condition,
  levels = c("Control", "ZIKV")
)

phenotype <- phenotype %>%
  dplyr::mutate(
    Condition_Label = dplyr::case_when(
      Condition == "Control" ~ "Control",
      Condition == "ZIKV" ~ "ZIKV infected"
    ),
    Condition_Label = factor(
      Condition_Label,
      levels = c("Control", "ZIKV infected")
    ),
    Plot_Label = paste0(SampleID, "\n", Condition_Label)
  )

rownames(phenotype) <- phenotype$geo_accession

write.csv(
  phenotype,
  "Results/Annotation/GSE146423_Final_6_Sample_Phenotype.csv",
  row.names = FALSE
)

cat("Final phenotype table:\n")
print(phenotype)

cat("\nPhenotype summary:\n")
print(table(phenotype$Condition, useNA = "ifany"))

if (nrow(phenotype) != 6) {
  stop("Expected 6 samples, but detected: ", nrow(phenotype))
}

if (any(is.na(phenotype$Condition))) {
  stop("Some samples could not be assigned to Control/ZIKV.")
}

cat("\n✓ Phenotype correctly defined for GSE146423.\n\n")


# ==============================================================================
# 4. COLORS FOR PUBLICATION FIGURES
# ==============================================================================

condition_colors <- c(
  "Control" = "#4DBBD5",
  "ZIKV infected" = "#E64B35"
)

deg_colors <- c(
  "Up in ZIKV" = "#D73027",
  "Down in ZIKV" = "#4575B4",
  "Not significant" = "grey75"
)

heatmap_annotation_colors <- list(
  Condition_Label = condition_colors,
  Cell_Line = c("A549" = "#7E6148")
)


# ==============================================================================
# 5. DOWNLOAD NCBI-GENERATED RAW COUNTS AND ANNOTATION
# ==============================================================================

cat("[4] Downloading NCBI/GEO-generated raw counts and annotation...\n")

options(timeout = 300000)

raw_counts_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE146423&file=GSE146423_raw_counts_GRCh38.p13_NCBI.tsv.gz&format=file&type=rnaseq_counts"

annotation_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?file=Human.GRCh38.p13.annot.tsv.gz&format=file&type=rnaseq_counts"

raw_counts_file <- "Results/Raw/GSE146423_raw_counts_GRCh38.p13_NCBI.tsv.gz"
annotation_file <- "Results/Raw/Human.GRCh38.p13.annot.tsv.gz"

if (!file.exists(raw_counts_file)) {
  download.file(
    url = raw_counts_url,
    destfile = raw_counts_file,
    mode = "wb"
  )
} else {
  cat("Using existing raw count file:", raw_counts_file, "\n")
}

if (!file.exists(annotation_file)) {
  download.file(
    url = annotation_url,
    destfile = annotation_file,
    mode = "wb"
  )
} else {
  cat("Using existing annotation file:", annotation_file, "\n")
}

if (!file.exists(raw_counts_file)) {
  stop("Raw count file was not downloaded/found.")
}

if (!file.exists(annotation_file)) {
  stop("Annotation file was not downloaded/found.")
}

cat("✓ Raw counts and annotation downloaded/found.\n\n")


# ==============================================================================
# 6. READ RAW COUNT MATRIX
# ==============================================================================

cat("[5] Reading raw count matrix...\n")

counts_raw <- data.table::fread(
  raw_counts_file,
  data.table = FALSE
)

write.csv(
  data.frame(Column_Names = colnames(counts_raw)),
  "Results/Annotation/GSE146423_Raw_Counts_Column_Names.csv",
  row.names = FALSE
)

cat("Raw count table dimensions:\n")
cat("Rows:", nrow(counts_raw), "\n")
cat("Columns:", ncol(counts_raw), "\n\n")
cat("Column names:\n")
print(colnames(counts_raw))
cat("\n")

gene_col_candidates <- c(
  "GeneID", "GeneId", "gene_id", "gene", "Gene", "ID",
  "EntrezID", "ENTREZID", "geneID", "geneid"
)

gene_col <- intersect(gene_col_candidates, colnames(counts_raw))

if (length(gene_col) == 0) {
  gene_col <- colnames(counts_raw)[1]
  warning("Could not identify gene ID column by name. Using first column: ", gene_col)
} else {
  gene_col <- gene_col[1]
}

cat("Selected gene ID column:", gene_col, "\n")

expected_samples <- phenotype$geo_accession

missing_samples <- setdiff(expected_samples, colnames(counts_raw))

if (length(missing_samples) > 0) {
  cat("Available columns in raw count table:\n")
  print(colnames(counts_raw))
  stop("Missing expected sample columns: ", paste(missing_samples, collapse = ", "))
}

count_df <- counts_raw[, c(gene_col, expected_samples)]
colnames(count_df)[1] <- "gene_id"

count_df$gene_id <- as.character(count_df$gene_id)

for (cc in expected_samples) {
  count_df[[cc]] <- suppressWarnings(as.numeric(count_df[[cc]]))
}

count_df <- count_df[!is.na(count_df$gene_id), ]
count_df <- count_df[count_df$gene_id != "", ]

count_df <- count_df %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(expected_samples), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

count_mat <- as.data.frame(count_df)
rownames(count_mat) <- count_mat$gene_id
count_mat$gene_id <- NULL

count_mat <- as.matrix(count_mat)
mode(count_mat) <- "numeric"
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)
mode(count_mat) <- "integer"

count_mat <- count_mat[, expected_samples, drop = FALSE]

if (!all(colnames(count_mat) == rownames(phenotype))) {
  cat("Count matrix columns:\n")
  print(colnames(count_mat))
  
  cat("Phenotype rownames:\n")
  print(rownames(phenotype))
  
  stop("Sample order mismatch between count matrix and phenotype.")
}

write.csv(
  count_mat,
  "Results/Counts/GSE146423_RAW_NCBI_Count_Matrix.csv"
)

cat("Raw count matrix dimensions:\n")
cat("Genes:", nrow(count_mat), "\n")
cat("Samples:", ncol(count_mat), "\n\n")
cat("Library sizes:\n")
print(colSums(count_mat))

cat("\n✓ Raw count matrix built successfully.\n\n")


# ==============================================================================
# 7. READ GENE ANNOTATION
# ==============================================================================

cat("[6] Reading gene annotation...\n")

annot_raw <- data.table::fread(
  annotation_file,
  data.table = FALSE
)

write.csv(
  data.frame(Column_Names = colnames(annot_raw)),
  "Results/Annotation/GSE146423_Annotation_Column_Names.csv",
  row.names = FALSE
)

cat("Annotation table dimensions:\n")
cat("Rows:", nrow(annot_raw), "\n")
cat("Columns:", ncol(annot_raw), "\n\n")
cat("Annotation column names:\n")
print(colnames(annot_raw))
cat("\n")

annot_gene_col_candidates <- c(
  "GeneID", "GeneId", "gene_id", "gene", "Gene", "ID",
  "EntrezID", "ENTREZID", "geneID", "geneid"
)

annot_gene_col <- intersect(annot_gene_col_candidates, colnames(annot_raw))

if (length(annot_gene_col) == 0) {
  annot_gene_col <- colnames(annot_raw)[1]
  warning("Could not identify annotation gene ID column by name. Using first column: ", annot_gene_col)
} else {
  annot_gene_col <- annot_gene_col[1]
}

symbol_col_candidates <- c("Symbol", "symbol", "GeneSymbol", "gene_symbol", "SYMBOL", "Name", "name")
symbol_col <- intersect(symbol_col_candidates, colnames(annot_raw))

description_col_candidates <- c("Description", "description", "GeneName", "gene_name", "GENENAME", "FullName")
description_col <- intersect(description_col_candidates, colnames(annot_raw))

gene_type_col_candidates <- c("type", "Type", "GeneType", "gene_type", "gene_biotype", "biotype")
gene_type_col <- intersect(gene_type_col_candidates, colnames(annot_raw))

annot_df <- data.frame(
  gene_id = as.character(annot_raw[[annot_gene_col]]),
  stringsAsFactors = FALSE
)

if (length(symbol_col) > 0) {
  annot_df$SYMBOL <- as.character(annot_raw[[symbol_col[1]]])
} else {
  annot_df$SYMBOL <- annot_df$gene_id
}

if (length(description_col) > 0) {
  annot_df$GENENAME <- as.character(annot_raw[[description_col[1]]])
} else {
  annot_df$GENENAME <- NA_character_
}

if (length(gene_type_col) > 0) {
  annot_df$GeneType <- as.character(annot_raw[[gene_type_col[1]]])
} else {
  annot_df$GeneType <- NA_character_
}

annot_df <- annot_df %>%
  dplyr::filter(!is.na(gene_id), gene_id != "") %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)

write.csv(
  annot_df,
  "Results/Annotation/GSE146423_Gene_Annotation_Parsed.csv",
  row.names = FALSE
)

cat("✓ Gene annotation parsed.\n\n")


# ==============================================================================
# 8. RAW COUNT QC PLOTS
# ==============================================================================

cat("[7] Generating raw count QC plots...\n")

lib_size <- colSums(count_mat)

lib_df <- data.frame(
  geo_accession = names(lib_size),
  Library_Size = as.numeric(lib_size),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(phenotype, by = "geo_accession") %>%
  dplyr::mutate(
    Plot_Label = factor(Plot_Label, levels = phenotype$Plot_Label)
  )

p_lib <- ggplot(
  lib_df,
  aes(x = Plot_Label, y = Library_Size, fill = Condition_Label)
) +
  geom_col(width = 0.72, color = "black", linewidth = 0.35) +
  geom_text(
    aes(label = scales::comma(round(Library_Size / 1e6, 1), suffix = "M")),
    vjust = -0.35,
    size = 4.2,
    fontface = "bold"
  ) +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(
    labels = scales::label_number(scale = 1e-6, suffix = "M"),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Library size per sample",
    subtitle = "NCBI/GEO-generated raw gene counts",
    x = "Sample",
    y = "Total reads assigned to genes",
    fill = "Group"
  ) +
  pub_theme

save_pub_plot(
  p_lib,
  "Results/QC/01_Publication_Library_Size",
  width = 8.5,
  height = 6
)

raw_long <- as.data.frame(log2(count_mat + 1)) %>%
  tibble::rownames_to_column("gene_id") %>%
  tidyr::pivot_longer(
    cols = -gene_id,
    names_to = "geo_accession",
    values_to = "log2_count"
  ) %>%
  dplyr::left_join(phenotype, by = "geo_accession") %>%
  dplyr::mutate(
    Plot_Label = factor(Plot_Label, levels = phenotype$Plot_Label)
  )

p_box_raw <- ggplot(
  raw_long,
  aes(x = Plot_Label, y = log2_count, fill = Condition_Label)
) +
  geom_boxplot(
    outlier.size = 0.15,
    linewidth = 0.45,
    width = 0.7,
    color = "black"
  ) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Raw count distribution",
    subtitle = expression(log[2]~"(raw counts + 1)"),
    x = "Sample",
    y = expression(log[2]~"raw count"),
    fill = "Group"
  ) +
  pub_theme

save_pub_plot(
  p_box_raw,
  "Results/QC/02_Publication_Raw_Count_Boxplot",
  width = 8.5,
  height = 6
)

cat("✓ Raw count QC plots completed.\n\n")


# ==============================================================================
# 9. LOW-EXPRESSION FILTERING
# ==============================================================================

cat("[8] Filtering low-expression genes...\n")

keep <- rowSums(count_mat >= MIN_COUNT) >= MIN_SAMPLES

filtered_count_mat <- count_mat[keep, , drop = FALSE]

write.csv(
  filtered_count_mat,
  "Results/Counts/GSE146423_Filtered_NCBI_Raw_Count_Matrix.csv"
)

filter_summary <- data.frame(
  Dataset = GSE_ID,
  Total_genes_before_filtering = nrow(count_mat),
  Total_genes_after_filtering = nrow(filtered_count_mat),
  Removed_genes = nrow(count_mat) - nrow(filtered_count_mat),
  Min_count = MIN_COUNT,
  Min_samples = MIN_SAMPLES
)

write.csv(
  filter_summary,
  "Results/Tables/GSE146423_Filtering_Summary.csv",
  row.names = FALSE
)

print(filter_summary)
cat("\n✓ Filtering completed.\n\n")


# ==============================================================================
# 10. DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS
# ==============================================================================

cat("[9] Running DESeq2 analysis...\n")

dds <- DESeqDataSetFromMatrix(
  countData = filtered_count_mat,
  colData = phenotype,
  design = ~ Condition
)

dds$Condition <- relevel(dds$Condition, ref = "Control")

dds <- DESeq(dds)

cat("DESeq2 result names:\n")
print(resultsNames(dds))

norm_counts <- counts(dds, normalized = TRUE)

write.csv(
  norm_counts,
  "Results/Counts/GSE146423_DESeq2_Normalized_Counts.csv"
)

vsd <- vst(dds, blind = FALSE)

vst_mat <- assay(vsd)

write.csv(
  vst_mat,
  "Results/Counts/GSE146423_VST_Counts.csv"
)

cat("✓ DESeq2 completed successfully.\n\n")


# ==============================================================================
# 11. NORMALIZATION QC PLOTS
# ==============================================================================

cat("[10] Generating normalization QC plots...\n")

norm_long <- as.data.frame(log2(norm_counts + 1)) %>%
  tibble::rownames_to_column("gene_id") %>%
  tidyr::pivot_longer(
    cols = -gene_id,
    names_to = "geo_accession",
    values_to = "log2_normalized_count"
  ) %>%
  dplyr::left_join(phenotype, by = "geo_accession") %>%
  dplyr::mutate(
    Plot_Label = factor(Plot_Label, levels = phenotype$Plot_Label)
  )

p_box_norm <- ggplot(
  norm_long,
  aes(x = Plot_Label, y = log2_normalized_count, fill = Condition_Label)
) +
  geom_boxplot(
    outlier.size = 0.15,
    linewidth = 0.45,
    width = 0.7,
    color = "black"
  ) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Normalized count distribution",
    subtitle = expression("DESeq2 normalized counts; " ~ log[2]~"(normalized counts + 1)"),
    x = "Sample",
    y = expression(log[2]~"normalized count"),
    fill = "Group"
  ) +
  pub_theme

save_pub_plot(
  p_box_norm,
  "Results/QC/03_Publication_DESeq2_Normalized_Count_Boxplot",
  width = 8.5,
  height = 6
)

# ------------------------------------------------------------------------------
# PCA plot - duplicate-column safe
# ------------------------------------------------------------------------------

pca_data <- DESeq2::plotPCA(
  vsd,
  intgroup = "Condition",
  returnData = TRUE
)

percentVar <- round(100 * attr(pca_data, "percentVar"))

if ("name" %in% colnames(pca_data)) {
  pca_data$geo_accession <- as.character(pca_data$name)
} else if (!"geo_accession" %in% colnames(pca_data)) {
  pca_data$geo_accession <- rownames(pca_data)
}

pca_data <- pca_data %>%
  dplyr::select(
    PC1,
    PC2,
    geo_accession
  ) %>%
  dplyr::left_join(phenotype, by = "geo_accession")

p_pca <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = Condition_Label, label = SampleID)
) +
  geom_point(size = 5.2, alpha = 0.95) +
  ggrepel::geom_text_repel(
    size = 5,
    fontface = "bold",
    box.padding = 0.55,
    point.padding = 0.35,
    segment.color = "grey35",
    segment.size = 0.35,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "Principal component analysis",
    subtitle = "A549 cells: ZIKV infected vs Control",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Group"
  ) +
  pub_theme +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.position = "top",
    legend.direction = "horizontal"
  )

save_pub_plot(
  p_pca,
  "Results/QC/04_Publication_PCA_ZIKV_vs_Control",
  width = 7.5,
  height = 6.5
)

sample_dists <- dist(t(vst_mat))
sample_dist_mat <- as.matrix(sample_dists)

rownames(sample_dist_mat) <- phenotype$Plot_Label[
  match(rownames(sample_dist_mat), phenotype$geo_accession)
]

colnames(sample_dist_mat) <- phenotype$Plot_Label[
  match(colnames(sample_dist_mat), phenotype$geo_accession)
]

ann_col <- phenotype[, c("Condition_Label", "Cell_Line")]
rownames(ann_col) <- phenotype$Plot_Label

save_pheatmap_pub(
  plot_function = function() {
    pheatmap(
      sample_dist_mat,
      annotation_col = ann_col,
      annotation_row = ann_col,
      annotation_colors = heatmap_annotation_colors,
      color = colorRampPalette(rev(brewer.pal(n = 9, name = "YlGnBu")))(100),
      border_color = NA,
      fontsize = 13,
      fontsize_row = 11,
      fontsize_col = 11,
      angle_col = 45,
      main = "Sample-to-sample distance"
    )
  },
  filename_base = "Results/QC/05_Publication_Sample_Distance_Heatmap",
  width = 8,
  height = 7.5
)

sample_cor_mat <- cor(vst_mat, method = "pearson")

rownames(sample_cor_mat) <- phenotype$Plot_Label[
  match(rownames(sample_cor_mat), phenotype$geo_accession)
]

colnames(sample_cor_mat) <- phenotype$Plot_Label[
  match(colnames(sample_cor_mat), phenotype$geo_accession)
]

save_pheatmap_pub(
  plot_function = function() {
    pheatmap(
      sample_cor_mat,
      annotation_col = ann_col,
      annotation_row = ann_col,
      annotation_colors = heatmap_annotation_colors,
      color = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(100),
      border_color = NA,
      fontsize = 13,
      fontsize_row = 11,
      fontsize_col = 11,
      angle_col = 45,
      display_numbers = TRUE,
      number_format = "%.2f",
      main = "Sample correlation heatmap"
    )
  },
  filename_base = "Results/QC/06_Publication_Sample_Correlation_Heatmap",
  width = 8,
  height = 7.5
)

save_base_pub(
  plot_function = function() {
    plotDispEsts(
      dds,
      main = "DESeq2 dispersion estimates"
    )
  },
  filename_base = "Results/QC/07_Publication_DESeq2_Dispersion_Plot",
  width = 7,
  height = 6
)

cat("✓ Normalization QC plots completed.\n\n")


# ==============================================================================
# 12. EXTRACT AND ANNOTATE DEG RESULTS
# ==============================================================================

cat("[11] Extracting and annotating DEG results...\n")

res <- results(
  dds,
  contrast = c("Condition", "ZIKV", "Control"),
  alpha = PADJ_CUTOFF
)

coef_name <- grep(
  "Condition.*ZIKV.*Control|ZIKV.*vs.*Control",
  resultsNames(dds),
  value = TRUE
)

if (length(coef_name) == 0) {
  coef_name <- resultsNames(dds)[length(resultsNames(dds))]
}

cat("Using coefficient for LFC shrinkage:", coef_name[1], "\n")

res_shrunk <- tryCatch(
  {
    lfcShrink(
      dds,
      coef = coef_name[1],
      type = "apeglm"
    )
  },
  error = function(e) {
    cat("apeglm shrinkage failed. Using unshrunken DESeq2 results.\n")
    cat("Reason:", e$message, "\n")
    res
  }
)

res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("gene_id")

res_unshrunk_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::select(gene_id, stat, pvalue, padj)

res_df <- res_df %>%
  dplyr::select(-dplyr::any_of(c("stat", "pvalue", "padj"))) %>%
  dplyr::left_join(res_unshrunk_df, by = "gene_id") %>%
  dplyr::arrange(padj)

res_df$gene_id <- as.character(res_df$gene_id)

res_df <- res_df %>%
  dplyr::left_join(annot_df, by = "gene_id")

missing_symbol <- is.na(res_df$SYMBOL) | res_df$SYMBOL == ""

if (sum(missing_symbol, na.rm = TRUE) > 0) {
  entrez_keys <- unique(res_df$gene_id[missing_symbol])
  entrez_keys <- entrez_keys[!is.na(entrez_keys) & entrez_keys != ""]
  
  if (length(entrez_keys) > 0) {
    org_map <- tryCatch(
      {
        AnnotationDbi::select(
          org.Hs.eg.db,
          keys = entrez_keys,
          keytype = "ENTREZID",
          columns = c("SYMBOL", "GENENAME")
        ) %>%
          dplyr::rename(gene_id = ENTREZID) %>%
          dplyr::distinct(gene_id, .keep_all = TRUE)
      },
      error = function(e) {
        NULL
      }
    )
    
    if (!is.null(org_map)) {
      res_df <- res_df %>%
        dplyr::left_join(org_map, by = "gene_id", suffix = c("", ".org")) %>%
        dplyr::mutate(
          SYMBOL = ifelse(is.na(SYMBOL) | SYMBOL == "", SYMBOL.org, SYMBOL),
          GENENAME = ifelse(is.na(GENENAME) | GENENAME == "", GENENAME.org, GENENAME)
        ) %>%
        dplyr::select(-dplyr::any_of(c("SYMBOL.org", "GENENAME.org")))
    }
  }
}

res_df <- res_df %>%
  dplyr::mutate(
    SYMBOL = ifelse(is.na(SYMBOL) | SYMBOL == "", gene_id, SYMBOL),
    ENTREZID = gene_id,
    Regulation = dplyr::case_when(
      !is.na(padj) & padj < PADJ_CUTOFF & log2FoldChange >= LFC_CUTOFF  ~ "Up in ZIKV",
      !is.na(padj) & padj < PADJ_CUTOFF & log2FoldChange <= -LFC_CUTOFF ~ "Down in ZIKV",
      TRUE ~ "Not significant"
    ),
    Regulation = factor(
      Regulation,
      levels = c("Up in ZIKV", "Down in ZIKV", "Not significant")
    ),
    minusLog10P = -log10(pvalue),
    minusLog10FDR = -log10(padj)
  )

res_df$minusLog10P[is.infinite(res_df$minusLog10P)] <- max(
  res_df$minusLog10P[is.finite(res_df$minusLog10P)],
  na.rm = TRUE
)

res_df$minusLog10FDR[is.infinite(res_df$minusLog10FDR)] <- max(
  res_df$minusLog10FDR[is.finite(res_df$minusLog10FDR)],
  na.rm = TRUE
)

write.csv(
  res_df,
  "Results/Tables/GSE146423_DESeq2_All_Genes_ZIKV_vs_Control.csv",
  row.names = FALSE
)

sig_deg <- res_df %>%
  dplyr::filter(
    !is.na(padj),
    padj < PADJ_CUTOFF,
    abs(log2FoldChange) >= LFC_CUTOFF
  )

deg_up <- sig_deg %>%
  dplyr::filter(log2FoldChange >= LFC_CUTOFF)

deg_down <- sig_deg %>%
  dplyr::filter(log2FoldChange <= -LFC_CUTOFF)

write.csv(
  sig_deg,
  "Results/Tables/GSE146423_DEG_FDR_0.05_log2FC_1.csv",
  row.names = FALSE
)

write.csv(
  deg_up,
  "Results/Tables/GSE146423_DEG_UP_in_ZIKV_FDR_0.05_log2FC_1.csv",
  row.names = FALSE
)

write.csv(
  deg_down,
  "Results/Tables/GSE146423_DEG_DOWN_in_ZIKV_FDR_0.05_log2FC_1.csv",
  row.names = FALSE
)

summary_table <- data.frame(
  Dataset = GSE_ID,
  Contrast = "ZIKV_vs_Control",
  Reference_group = "Control",
  Test_group = "ZIKV",
  Count_type = "NCBI/GEO-generated raw RNA-seq gene counts",
  Total_genes_raw = nrow(count_mat),
  Total_genes_after_filtering = nrow(filtered_count_mat),
  FDR_cutoff = PADJ_CUTOFF,
  log2FC_cutoff = LFC_CUTOFF,
  Significant_DEG = nrow(sig_deg),
  Up_in_ZIKV = nrow(deg_up),
  Down_in_ZIKV = nrow(deg_down),
  Log2FC_interpretation = "Positive log2FC means higher expression in ZIKV infected cells"
)

write.csv(
  summary_table,
  "Results/Tables/GSE146423_DEG_Summary.csv",
  row.names = FALSE
)

cat("DEG summary:\n")
print(summary_table)
cat("\n✓ DEG extraction and annotation completed.\n\n")


# ==============================================================================
# 13. PUBLICATION-READY DEG PLOTS
# ==============================================================================

cat("[12] Generating publication-ready DEG plots...\n")

deg_count_df <- res_df %>%
  dplyr::count(Regulation) %>%
  dplyr::mutate(
    Regulation = factor(Regulation, levels = levels(res_df$Regulation))
  )

p_deg_count <- ggplot(
  deg_count_df,
  aes(x = Regulation, y = n, fill = Regulation)
) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.4, fontface = "bold", size = 4) +
  scale_fill_manual(values = deg_colors) +
  labs(
    title = "Differentially expressed genes",
    subtitle = "FDR < 0.05 and |log2FC| >= 1",
    x = NULL,
    y = "Number of genes"
  ) +
  pub_theme +
  theme(legend.position = "none")

save_pub_plot(
  p_deg_count,
  "Results/Plots/01_Publication_DEG_Count_Barplot",
  width = 7,
  height = 6
)

save_base_pub(
  plot_function = function() {
    plotMA(
      res_shrunk,
      ylim = c(-6, 6),
      main = "MA plot: ZIKV infected vs Control"
    )
    abline(
      h = c(-LFC_CUTOFF, LFC_CUTOFF),
      col = "dodgerblue4",
      lty = 2,
      lwd = 1.5
    )
  },
  filename_base = "Results/Plots/02_Publication_MA_Plot_ZIKV_vs_Control",
  width = 7,
  height = 6
)

# Data-driven labels only: no manually selected showcase genes.
top_up_labels <- res_df %>%
  dplyr::filter(Regulation == "Up in ZIKV") %>%
  dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::slice_head(n = 12)

top_down_labels <- res_df %>%
  dplyr::filter(Regulation == "Down in ZIKV") %>%
  dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::slice_head(n = 12)

label_df <- dplyr::bind_rows(top_up_labels, top_down_labels) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

write.csv(
  label_df,
  "Results/Tables/GSE146423_Genes_Labelled_In_Publication_Plots.csv",
  row.names = FALSE
)

volcano_counts <- res_df %>%
  dplyr::count(Regulation)

get_count <- function(group_name) {
  value <- volcano_counts$n[match(group_name, volcano_counts$Regulation)]
  ifelse(is.na(value), 0, value)
}

legend_labels <- c(
  "Up in ZIKV" = paste0("Up in ZIKV (n = ", get_count("Up in ZIKV"), ")"),
  "Down in ZIKV" = paste0("Down in ZIKV (n = ", get_count("Down in ZIKV"), ")"),
  "Not significant" = paste0("Not significant (n = ", get_count("Not significant"), ")")
)

p_volcano <- ggplot(
  res_df,
  aes(x = log2FoldChange, y = minusLog10FDR, color = Regulation)
) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_vline(
    xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = -log10(PADJ_CUTOFF),
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_color_manual(values = deg_colors, labels = legend_labels, drop = FALSE) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = SYMBOL),
    size = 3.5,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey35",
    segment.size = 0.35
  ) +
  labs(
    title = "Volcano plot: ZIKV vs Control",
    subtitle = "A549 cells; DESeq2 raw-count analysis; data-driven labels",
    x = expression(log[2]~fold~change~"(ZIKV / Control)"),
    y = expression(-log[10]~FDR),
    color = "Regulation"
  ) +
  pub_theme +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.position = "right"
  )

save_pub_plot(
  p_volcano,
  "Results/Plots/03_Publication_Volcano_ZIKV_vs_Control",
  width = 10,
  height = 7.5
)

make_safe_pdf(
  "Results/Plots/04_Publication_EnhancedVolcano_ZIKV_vs_Control.pdf",
  {
    print(
      EnhancedVolcano(
        res_df,
        lab = res_df$SYMBOL,
        x = "log2FoldChange",
        y = "padj",
        selectLab = label_df$SYMBOL,
        title = "GSE146423: ZIKV vs Control",
        subtitle = "A549 cells; NCBI raw counts analyzed with DESeq2",
        caption = "Cutoff: FDR < 0.05 and |log2FC| >= 1",
        pCutoff = PADJ_CUTOFF,
        FCcutoff = LFC_CUTOFF,
        pointSize = 2.4,
        labSize = 4.0,
        boxedLabels = TRUE,
        colAlpha = 0.85,
        legendPosition = "right",
        drawConnectors = TRUE,
        widthConnectors = 0.5,
        colConnectors = "grey30"
      )
    )
  },
  width = 10,
  height = 8
)

png(
  "Results/Plots/04_Publication_EnhancedVolcano_ZIKV_vs_Control.png",
  width = 10,
  height = 8,
  units = "in",
  res = 600
)

print(
  EnhancedVolcano(
    res_df,
    lab = res_df$SYMBOL,
    x = "log2FoldChange",
    y = "padj",
    selectLab = label_df$SYMBOL,
    title = "GSE146423: ZIKV vs Control",
    subtitle = "A549 cells; NCBI raw counts analyzed with DESeq2",
    caption = "Cutoff: FDR < 0.05 and |log2FC| >= 1",
    pCutoff = PADJ_CUTOFF,
    FCcutoff = LFC_CUTOFF,
    pointSize = 2.4,
    labSize = 4.0,
    boxedLabels = TRUE,
    colAlpha = 0.85,
    legendPosition = "right",
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = "grey30"
  )
)

dev.off()

top_up <- res_df %>%
  dplyr::filter(Regulation == "Up in ZIKV") %>%
  dplyr::arrange(desc(log2FoldChange)) %>%
  dplyr::slice_head(n = 15)

top_down <- res_df %>%
  dplyr::filter(Regulation == "Down in ZIKV") %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::slice_head(n = 15)

top_bar_df <- dplyr::bind_rows(top_down, top_up) %>%
  dplyr::mutate(
    Gene_Label = ifelse(is.na(SYMBOL), gene_id, SYMBOL),
    Gene_Label = factor(Gene_Label, levels = Gene_Label[order(log2FoldChange)])
  )

write.csv(
  top_bar_df,
  "Results/Tables/GSE146423_Top_Up_Down_Genes_For_Barplot.csv",
  row.names = FALSE
)

if (nrow(top_bar_df) >= 1) {
  p_top_bar <- ggplot(
    top_bar_df,
    aes(x = Gene_Label, y = log2FoldChange, fill = Regulation)
  ) +
    geom_col(color = "black", linewidth = 0.25, width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = deg_colors) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(
      title = "Top differentially expressed genes",
      subtitle = "Top genes ranked by log2 fold change",
      x = NULL,
      y = "log2 fold change",
      fill = "Regulation"
    ) +
    pub_theme
  
  save_pub_plot(
    p_top_bar,
    "Results/Plots/05_Publication_Top_Up_Down_DEG_Barplot",
    width = 8,
    height = 9
  )
}

top_lollipop_df <- res_df %>%
  dplyr::filter(!is.na(padj), padj < PADJ_CUTOFF) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::mutate(
    Gene_Label = ifelse(is.na(SYMBOL), gene_id, SYMBOL),
    Gene_Label = factor(Gene_Label, levels = rev(Gene_Label))
  )

write.csv(
  top_lollipop_df,
  "Results/Tables/GSE146423_Top30_DEG_For_Lollipop.csv",
  row.names = FALSE
)

if (nrow(top_lollipop_df) >= 1) {
  p_lollipop <- ggplot(
    top_lollipop_df,
    aes(x = minusLog10FDR, y = Gene_Label, color = Regulation)
  ) +
    geom_segment(
      aes(x = 0, xend = minusLog10FDR, y = Gene_Label, yend = Gene_Label),
      linewidth = 0.7,
      color = "grey60"
    ) +
    geom_point(size = 3.2) +
    scale_color_manual(values = deg_colors) +
    labs(
      title = "Top 30 genes ranked by FDR",
      subtitle = "Higher value indicates stronger statistical significance",
      x = "-log10 FDR",
      y = NULL,
      color = "Regulation"
    ) +
    pub_theme
  
  save_pub_plot(
    p_lollipop,
    "Results/Plots/06_Publication_Top30_DEG_Lollipop",
    width = 8,
    height = 9
  )
}

rank_df <- res_df %>%
  dplyr::arrange(desc(log2FoldChange)) %>%
  dplyr::mutate(Rank = dplyr::row_number())

p_rank <- ggplot(
  rank_df,
  aes(x = Rank, y = log2FoldChange, color = Regulation)
) +
  geom_point(alpha = 0.75, size = 1.4) +
  geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed")) +
  scale_color_manual(values = deg_colors) +
  labs(
    title = "Ranked log2 fold-change profile",
    subtitle = "Genes ordered from highest to lowest log2FC",
    x = "Gene rank",
    y = "log2 fold change",
    color = "Regulation"
  ) +
  pub_theme +
  theme(axis.text.x = element_text(angle = 0))

save_pub_plot(
  p_rank,
  "Results/Plots/07_Publication_Ranked_log2FC_Profile",
  width = 8,
  height = 6
)

cat("✓ Publication-ready DEG plots completed.\n\n")


# ==============================================================================
# 14. HEATMAPS OF TOP DEGS
# ==============================================================================

cat("[13] Generating DEG heatmaps...\n")

make_deg_heatmap <- function(n_genes, filename_base, width, height, fontsize_row) {
  
  if (nrow(sig_deg) < 2) {
    cat("Skipping heatmap: fewer than 2 significant DEGs.\n")
    return(NULL)
  }
  
  top_genes <- sig_deg %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = min(n_genes, nrow(sig_deg))) %>%
    dplyr::pull(gene_id)
  
  top_genes <- intersect(top_genes, rownames(vst_mat))
  
  if (length(top_genes) < 2) {
    cat("Skipping heatmap: fewer than 2 valid genes.\n")
    return(NULL)
  }
  
  heat_mat <- vst_mat[top_genes, , drop = FALSE]
  heat_mat_z <- t(scale(t(heat_mat)))
  heat_mat_z[is.na(heat_mat_z)] <- 0
  
  symbol_lookup <- res_df %>%
    dplyr::select(gene_id, SYMBOL) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)
  
  new_rownames <- symbol_lookup$SYMBOL[
    match(rownames(heat_mat_z), symbol_lookup$gene_id)
  ]
  
  rownames(heat_mat_z) <- make.unique(
    ifelse(is.na(new_rownames), rownames(heat_mat_z), new_rownames)
  )
  
  colnames(heat_mat_z) <- phenotype$Plot_Label[
    match(colnames(heat_mat_z), phenotype$geo_accession)
  ]
  
  ann_col_heat <- phenotype[, c("Condition_Label", "Cell_Line")]
  rownames(ann_col_heat) <- phenotype$Plot_Label
  
  save_pheatmap_pub(
    plot_function = function() {
      pheatmap(
        heat_mat_z,
        annotation_col = ann_col_heat,
        annotation_colors = heatmap_annotation_colors,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize = 12,
        fontsize_row = fontsize_row,
        fontsize_col = 10,
        angle_col = 45,
        border_color = NA,
        color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
        breaks = seq(-2.5, 2.5, length.out = 101),
        main = paste0("Top ", length(top_genes), " differentially expressed genes")
      )
    },
    filename_base = filename_base,
    width = width,
    height = height
  )
}

make_deg_heatmap(
  n_genes = TOP_N_MAIN_HEATMAP,
  filename_base = "Results/Plots/08_Publication_Top30_DEG_Heatmap",
  width = 8.5,
  height = 10,
  fontsize_row = 10
)

make_deg_heatmap(
  n_genes = TOP_N_SUPP_HEATMAP,
  filename_base = "Results/Plots/09_Supplementary_Top50_DEG_Heatmap",
  width = 8.5,
  height = 12,
  fontsize_row = 8
)

cat("✓ DEG heatmaps completed.\n\n")


# ==============================================================================
# 15. FUNCTIONAL ENRICHMENT
# ==============================================================================

cat("[14] Running GO/KEGG enrichment...\n")

main_deg <- sig_deg

up_deg <- main_deg %>%
  dplyr::filter(log2FoldChange >= LFC_CUTOFF)

down_deg <- main_deg %>%
  dplyr::filter(log2FoldChange <= -LFC_CUTOFF)

write.csv(
  main_deg,
  "Results/Tables/GSE146423_Main_DEG_Used_For_Enrichment.csv",
  row.names = FALSE
)

run_enrichment <- function(deg_table, prefix_name) {
  
  gene_table <- deg_table %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
  
  if (nrow(gene_table) < 5) {
    cat("Skipping enrichment for", prefix_name, ": fewer than 5 Entrez genes.\n")
    return(NULL)
  }
  
  entrez_ids <- unique(gene_table$ENTREZID)
  
  cat("Running enrichment for:", prefix_name, "\n")
  cat("Entrez genes:", length(entrez_ids), "\n")
  
  ego_bp <- tryCatch(
    enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) {
      cat("GO BP failed for", prefix_name, ":", e$message, "\n")
      NULL
    }
  )
  
  ego_mf <- tryCatch(
    enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) {
      cat("GO MF failed for", prefix_name, ":", e$message, "\n")
      NULL
    }
  )
  
  ego_cc <- tryCatch(
    enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) {
      cat("GO CC failed for", prefix_name, ":", e$message, "\n")
      NULL
    }
  )
  
  ekegg <- tryCatch(
    enrichKEGG(
      gene = entrez_ids,
      organism = "hsa",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    ),
    error = function(e) {
      cat("KEGG failed for", prefix_name, ":", e$message, "\n")
      NULL
    }
  )
  
  if (!is.null(ekegg)) {
    ekegg <- tryCatch(
      setReadable(
        ekegg,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID"
      ),
      error = function(e) ekegg
    )
  }
  
  enrichment_list <- list(
    GO_BP = ego_bp,
    GO_MF = ego_mf,
    GO_CC = ego_cc,
    KEGG = ekegg
  )
  
  for (nm in names(enrichment_list)) {
    
    obj <- enrichment_list[[nm]]
    
    if (is.null(obj)) {
      next
    }
    
    res_enrich <- as.data.frame(obj)
    
    write.csv(
      res_enrich,
      paste0("Results/Enrichment/", prefix_name, "_", nm, "_Enrichment.csv"),
      row.names = FALSE
    )
    
    if (nrow(res_enrich) > 0) {
      publication_enrich_dotplot(
        obj,
        paste0(prefix_name, " ", nm, " enrichment"),
        paste0("Results/Enrichment/", prefix_name, "_", nm, "_Publication_Dotplot"),
        show_n = 15
      )
    }
  }
  
  return(enrichment_list)
}

enrich_all <- run_enrichment(
  deg_table = main_deg,
  prefix_name = "GSE146423_All_Significant_DEG"
)

enrich_up <- run_enrichment(
  deg_table = up_deg,
  prefix_name = "GSE146423_Up_in_ZIKV_DEG"
)

enrich_down <- run_enrichment(
  deg_table = down_deg,
  prefix_name = "GSE146423_Down_in_ZIKV_DEG"
)

cat("✓ GO/KEGG enrichment completed.\n\n")


# ==============================================================================
# 16. GSEA ANALYSIS
# ==============================================================================

cat("[15] Running GSEA GO/KEGG...\n")

rank_df <- res_df %>%
  dplyr::filter(!is.na(ENTREZID), !is.na(stat)) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_max(
    order_by = abs(stat),
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

gene_rank <- rank_df$stat
names(gene_rank) <- rank_df$ENTREZID
gene_rank <- sort(gene_rank, decreasing = TRUE)

write.csv(
  data.frame(
    ENTREZID = names(gene_rank),
    statistic = as.numeric(gene_rank)
  ),
  "Results/Enrichment/GSE146423_GSEA_Ranked_Gene_List.csv",
  row.names = FALSE
)

if (length(gene_rank) >= 100) {
  
  gsea_go_bp <- tryCatch(
    {
      gseGO(
        geneList = gene_rank,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        eps = 0,
        verbose = FALSE
      )
    },
    error = function(e) {
      cat("GSEA GO BP failed:", e$message, "\n")
      NULL
    }
  )
  
  if (!is.null(gsea_go_bp)) {
    gsea_go_bp <- setReadable(
      gsea_go_bp,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID"
    )
    
    gsea_go_df_raw <- as.data.frame(gsea_go_bp)
    gsea_go_df_clean <- clean_gsea_df(gsea_go_df_raw)
    
    write.csv(
      gsea_go_df_raw,
      "Results/Enrichment/GSE146423_GSEA_GO_BP_ZIKV_vs_Control_RAW.csv",
      row.names = FALSE
    )
    
    write.csv(
      gsea_go_df_clean,
      "Results/Enrichment/GSE146423_GSEA_GO_BP_ZIKV_vs_Control_CLEAN.csv",
      row.names = FALSE
    )
    
    if (nrow(gsea_go_df_clean) > 0) {
      publication_enrich_dotplot(
        gsea_go_df_clean,
        "GSEA GO Biological Process: ZIKV vs Control",
        "Results/Enrichment/GSE146423_GSEA_GO_BP_Publication_Dotplot",
        show_n = 15
      )
      
      top_gsea_id <- gsea_go_df_clean$ID[1]
      
      gsea_top_plot <- tryCatch(
        gseaplot2(
          gsea_go_bp,
          geneSetID = top_gsea_id,
          title = gsea_go_df_clean$Description[1]
        ),
        error = function(e) {
          cat("Top GO GSEA plot failed:", e$message, "\n")
          NULL
        }
      )
      
      if (!is.null(gsea_top_plot)) {
        save_pub_plot(
          gsea_top_plot,
          "Results/Enrichment/GSE146423_GSEA_GO_BP_Top_Pathway",
          width = 9,
          height = 6.5
        )
      }
    }
  }
  
  gsea_kegg <- tryCatch(
    {
      gseKEGG(
        geneList = gene_rank,
        organism = "hsa",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        eps = 0,
        verbose = FALSE
      )
    },
    error = function(e) {
      cat("GSEA KEGG failed:", e$message, "\n")
      NULL
    }
  )
  
  if (!is.null(gsea_kegg)) {
    gsea_kegg_df_raw <- as.data.frame(gsea_kegg)
    gsea_kegg_df_clean <- clean_gsea_df(gsea_kegg_df_raw)
    
    write.csv(
      gsea_kegg_df_raw,
      "Results/Enrichment/GSE146423_GSEA_KEGG_ZIKV_vs_Control_RAW.csv",
      row.names = FALSE
    )
    
    write.csv(
      gsea_kegg_df_clean,
      "Results/Enrichment/GSE146423_GSEA_KEGG_ZIKV_vs_Control_CLEAN.csv",
      row.names = FALSE
    )
    
    if (nrow(gsea_kegg_df_clean) > 0) {
      publication_enrich_dotplot(
        gsea_kegg_df_clean,
        "GSEA KEGG: ZIKV vs Control",
        "Results/Enrichment/GSE146423_GSEA_KEGG_Publication_Dotplot",
        show_n = 15
      )
    }
  }
  
} else {
  cat("Skipping GSEA: fewer than 100 ranked genes.\n")
}

cat("✓ GSEA completed.\n\n")


# ==============================================================================
# 17. SAVE SESSION INFO
# ==============================================================================

sink("Results/Annotation/sessionInfo.txt")
print(sessionInfo())
sink()


# ==============================================================================
# 18. FINISH
# ==============================================================================

cat("\n============================================================\n")
cat("WORKFLOW COMPLETED SUCCESSFULLY\n")
cat("Dataset:", GSE_ID, "\n")
cat("Input type:\n")
cat("  NCBI/GEO-generated raw RNA-seq gene counts\n")
cat("Analysis type:\n")
cat("  Raw-count DESeq2 differential expression analysis\n")
cat("Comparison:\n")
cat("  ZIKV infected A549 cells vs Control A549 cells\n")
cat("Reference group:\n")
cat("  Control\n")
cat("Interpretation:\n")
cat("  Positive log2FoldChange = higher expression in ZIKV infected cells\n")
cat("\nMain output folders:\n")
cat("  Results/QC\n")
cat("  Results/Plots\n")
cat("  Results/Tables\n")
cat("  Results/Counts\n")
cat("  Results/Annotation\n")
cat("  Results/Enrichment\n")
cat("\nMain publication-ready plots:\n")
cat("  Results/QC/01_Publication_Library_Size.pdf/png\n")
cat("  Results/QC/02_Publication_Raw_Count_Boxplot.pdf/png\n")
cat("  Results/QC/03_Publication_DESeq2_Normalized_Count_Boxplot.pdf/png\n")
cat("  Results/QC/04_Publication_PCA_ZIKV_vs_Control.pdf/png\n")
cat("  Results/QC/05_Publication_Sample_Distance_Heatmap.pdf/png\n")
cat("  Results/QC/06_Publication_Sample_Correlation_Heatmap.pdf/png\n")
cat("  Results/QC/07_Publication_DESeq2_Dispersion_Plot.pdf/png\n")
cat("  Results/Plots/01_Publication_DEG_Count_Barplot.pdf/png\n")
cat("  Results/Plots/02_Publication_MA_Plot_ZIKV_vs_Control.pdf/png\n")
cat("  Results/Plots/03_Publication_Volcano_ZIKV_vs_Control.pdf/png\n")
cat("  Results/Plots/04_Publication_EnhancedVolcano_ZIKV_vs_Control.pdf/png\n")
cat("  Results/Plots/05_Publication_Top_Up_Down_DEG_Barplot.pdf/png\n")
cat("  Results/Plots/06_Publication_Top30_DEG_Lollipop.pdf/png\n")
cat("  Results/Plots/07_Publication_Ranked_log2FC_Profile.pdf/png\n")
cat("  Results/Plots/08_Publication_Top30_DEG_Heatmap.pdf/png\n")
cat("  Results/Plots/09_Supplementary_Top50_DEG_Heatmap.pdf/png\n")
cat("============================================================\n")

