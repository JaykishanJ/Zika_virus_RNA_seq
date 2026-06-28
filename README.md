# Zika Virus RNA-seq Meta-Analysis

This repository contains the code and methodology for a comprehensive RNA-seq meta-analysis investigating the transcriptional response of A549 cells to Zika Virus (ZIKV) infection. 

## Datasets
The analysis combines two independent datasets from NCBI GEO:
- **GSE146423**: RNA-seq of A549 cells (ZIKV vs Control), analyzed using NCBI-generated raw gene counts.
- **GSE265922**: RNA-seq of A549 cells (ZIKV vs Mock), analyzed using STAR unstranded reads.

## Project Structure
The project is organized as an RStudio project (`Zika_virus_wetlab.Rproj`) to ensure portability and reproducibility via the `here` package.

- `R/`: Contains shared utility functions (`utils.R`) and plotting functions (`plot_functions.R`) to keep the codebase DRY and maintain consistent publication-ready aesthetics.
- `GSE146423/`: Contains the specific differential expression analysis script for the GSE146423 dataset (`DEG_GSE146423.R`).
- `GSE265922/`: Contains the specific differential expression analysis script for the GSE265922 dataset (`DEG_GES265922.R`).
- `Meta_Analysis/`: Contains the script (`Meta_Analysis.R`) that merges the count matrices from both datasets and runs a unified DESeq2 model, controlling for batch effects between the two studies (`~ Dataset + Condition`).
- `setup_renv.R`: A script to initialize `renv` for dependency tracking.

## Reproducibility
This project uses `renv` to manage R package dependencies. 
1. Open `Zika_virus_wetlab.Rproj` in RStudio.
2. Run `source("setup_renv.R")` to initialize the environment and install necessary dependencies.

## Output Structure
Each analysis script (including the meta-analysis) generates its outputs locally within its respective directory under a `Results/` folder.
- `Results/Counts/`: Count matrices.
- `Results/QC/`: Quality control plots (PCA, dispersion, sample correlation).
- `Results/Tables/`: DESeq2 output tables, significant DEG lists, and GO/KEGG enrichment tables.
- `Results/Plots/`: Publication-ready visualizations (Volcano plots, MA plots, Heatmaps, etc.) in both high-res PNG and PDF formats.
