# ==============================================================================
# Plotting Functions (plot_functions.R)
# Contains shared plotting functions for consistent publication-ready aesthetics
# ==============================================================================

library(ggplot2)

pub_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16, color = "black"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6)
  )

make_safe_pdf <- function(filename, expr, width = 8, height = 6) {
  pdf(filename, width = width, height = height)
  
  tryCatch(
    expr,
    error = function(e) {
      plot.new()
      text(
        0.5,
        0.5,
        paste("Plot failed:\n", e$message),
        cex = 1
      )
      message("Plot failed: ", filename, " | ", e$message)
    },
    finally = dev.off()
  )
}

save_pub_plot <- function(plot_object, filename_base, width = 8, height = 6, dpi = 600) {
  print(plot_object)
  
  pdf_file <- paste0(filename_base, ".pdf")
  png_file <- paste0(filename_base, ".png")
  
  ggsave(
    filename = pdf_file,
    plot = plot_object,
    width = width,
    height = height,
    device = cairo_pdf
  )
  
  ggsave(
    filename = png_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

save_pheatmap_pub <- function(plot_function, filename_base, width = 8, height = 8, dpi = 600) {
  plot_function()
  
  pdf(
    file = paste0(filename_base, ".pdf"),
    width = width,
    height = height
  )
  plot_function()
  dev.off()
  
  png(
    filename = paste0(filename_base, ".png"),
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )
  plot_function()
  dev.off()
}

save_base_pub <- function(plot_function, filename_base, width = 8, height = 6, dpi = 600) {
  plot_function()
  
  pdf(
    file = paste0(filename_base, ".pdf"),
    width = width,
    height = height
  )
  plot_function()
  dev.off()
  
  png(
    filename = paste0(filename_base, ".png"),
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )
  plot_function()
  dev.off()
}

publication_enrich_dotplot <- function(enrich_object, title_text, filename_base, show_n = 15) {
  
  enrich_df <- as.data.frame(enrich_object)
  
  if (nrow(enrich_df) == 0) {
    cat("No enriched terms for:", title_text, "\n")
    return(NULL)
  }
  
  if ("Count" %in% colnames(enrich_df)) {
    enrich_df$Gene_Count_For_Plot <- enrich_df$Count
  } else if ("setSize" %in% colnames(enrich_df)) {
    enrich_df$Gene_Count_For_Plot <- enrich_df$setSize
  } else {
    enrich_df$Gene_Count_For_Plot <- seq_len(nrow(enrich_df))
  }
  
  if ("p.adjust" %in% colnames(enrich_df)) {
    enrich_df$Adjusted_P_For_Plot <- enrich_df$p.adjust
  } else if ("qvalues" %in% colnames(enrich_df)) {
    enrich_df$Adjusted_P_For_Plot <- enrich_df$qvalues
  } else if ("qvalue" %in% colnames(enrich_df)) {
    enrich_df$Adjusted_P_For_Plot <- enrich_df$qvalue
  } else if ("pvalue" %in% colnames(enrich_df)) {
    enrich_df$Adjusted_P_For_Plot <- enrich_df$pvalue
  } else {
    enrich_df$Adjusted_P_For_Plot <- seq_len(nrow(enrich_df))
  }
  
  enrich_df <- enrich_df %>%
    dplyr::arrange(Adjusted_P_For_Plot) %>%
    dplyr::slice_head(n = show_n) %>%
    dplyr::mutate(
      Description = stringr::str_wrap(Description, width = 42),
      Description = forcats::fct_reorder(Description, Gene_Count_For_Plot)
    )
  
  p <- ggplot(
    enrich_df,
    aes(
      x = Gene_Count_For_Plot,
      y = Description,
      size = Gene_Count_For_Plot,
      color = Adjusted_P_For_Plot
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_color_viridis_c(
      option = "plasma",
      direction = -1,
      name = "Adjusted p-value"
    ) +
    scale_size_continuous(
      name = "Gene count / set size",
      range = c(3, 8)
    ) +
    labs(
      title = title_text,
      x = "Gene count / gene-set size",
      y = NULL
    ) +
    pub_theme +
    theme(
      axis.text.x = element_text(angle = 0),
      axis.text.y = element_text(size = 11),
      legend.position = "right"
    )
  
  save_pub_plot(
    p,
    filename_base,
    width = 10.5,
    height = 7.5
  )
}
