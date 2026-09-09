#' Summary statistics for a Seurat object
#'
#' Computes cell count, total and per-cell UMI/gene statistics, and returns
#' both a formatted text label and a named list of values.
#'
#' @param seurat_obj A Seurat object.
#' @param title_prefix Optional string prepended to the summary label.
#' @return A list with `label` (character) and `stats` (named list of numeric
#'   values: n_cells, total_UMI, total_genes, avg_UMI_per_cell,
#'   avg_genes_per_cell, range_UMI, range_genes, avg_UMI_per_gene).
#' @export
generate_summary_stats <- function(seurat_obj, title_prefix = "") {
  meta <- seurat_obj@meta.data
  
  # Summary statistics
  n_cells <- ncol(seurat_obj)
  total_UMI <- sum(meta$nCount_RNA)
  total_genes <- sum(meta$nFeature_RNA)
  avg_UMI_per_cell <- round(mean(meta$nCount_RNA), 1)
  avg_genes_per_cell <- round(mean(meta$nFeature_RNA), 1)
  range_UMI <- range(meta$nCount_RNA)
  range_genes <- range(meta$nFeature_RNA)
  avg_UMI_per_gene <- round(total_UMI / total_genes, 2)
  
  # Create summary label
  summary_label <- paste(
    if(title_prefix != "") paste0(title_prefix, "\n") else "",
    "Cells:", n_cells,
    "\nUMIs: ", total_UMI,
    "\nGenes: ", total_genes,
    "\nAvg UMIs/cell: ", avg_UMI_per_cell,
    "\nRange UMIs/cell: ", paste(range_UMI, collapse = " - "),
    "\nAvg Genes/cell: ", avg_genes_per_cell,
    "\nRange Genes/cell: ", paste(range_genes, collapse = " - "),
    "\nAvg UMIs/gene: ", avg_UMI_per_gene,
    sep = ""
  )
  
  # Return both the label and individual stats
  return(list(
    label = summary_label,
    stats = list(
      n_cells = n_cells,
      total_UMI = total_UMI,
      total_genes = total_genes,
      avg_UMI_per_cell = avg_UMI_per_cell,
      avg_genes_per_cell = avg_genes_per_cell,
      range_UMI = range_UMI,
      range_genes = range_genes,
      avg_UMI_per_gene = avg_UMI_per_gene
    )
  ))
}

#' Stacked bar / alluvial plot of sample composition by cluster
#'
#' Shows the percentage each cluster contributes within each sample as a
#' stacked bar chart with alluvial flows connecting matching clusters across
#' samples.
#'
#' @param seurat_obj A Seurat object.
#' @param group_by Metadata column used for the fill grouping (default
#'   `"seurat_clusters"`).
#' @param sample_col Metadata column identifying samples.
#' @param sample_order Optional character vector specifying x-axis sample order.
#' @param title Plot title.
#' @param x_lab,y_lab Axis labels.
#' @param fill_lab Legend title for the fill; defaults to `group_by`.
#' @param show_percentages If `TRUE`, print percentage labels inside bars.
#' @param text_size Size of percentage labels.
#' @param bar_width Width of the bars.
#' @param flow_alpha Alpha for the alluvial flows.
#' @return A ggplot object.
#' @export
plot_sample_composition <- function(seurat_obj,
                                    group_by = "seurat_clusters",
                                    sample_col = "genotype",
                                    sample_order = NULL,
                                    title = "Sample Composition by Cluster",
                                    x_lab = "Sample",
                                    y_lab = "Percent Composition",
                                    fill_lab = NULL,
                                    show_percentages = TRUE,
                                    text_size = 3,
                                    bar_width = 0.5,
                                    flow_alpha = 0.5) {
  
  # Load required packages
  require(dplyr)
  require(ggplot2)
  require(ggalluvial)
  require(cowplot)
  
  # Extract metadata properly - convert to vectors
  cluster_data <- as.character(seurat_obj[[group_by]][, 1])
  sample_data <- as.character(seurat_obj[[sample_col]][, 1])
  
  # Create counts table
  cluster_counts <- as.data.frame(table(cluster_data, sample_data))
  colnames(cluster_counts) <- c("Cluster", "Sample", "Count")
  
  # Compute total cell count per sample
  sample_totals <- cluster_counts %>%
    group_by(Sample) %>%
    summarize(Total = sum(Count))
  
  # Calculate percentage composition
  cluster_percentages <- cluster_counts %>%
    left_join(sample_totals, by = "Sample") %>%
    mutate(Percent = (Count / Total) * 100)
  
  # Set sample order if provided
  if (!is.null(sample_order)) {
    cluster_percentages$Sample <- factor(cluster_percentages$Sample, levels = sample_order)
  }
  
  # Set fill label
  if (is.null(fill_lab)) {
    fill_lab <- group_by
  }
  
  # Create plot
  p <- ggplot(cluster_percentages, aes(x = Sample, y = Percent, fill = Cluster)) +
    geom_flow(aes(alluvium = Cluster), alpha = flow_alpha, curve_type = "linear", width = bar_width) +
    geom_col(width = bar_width) +
    cowplot::theme_minimal_hgrid() +
    labs(title = title, x = x_lab, y = y_lab, fill = fill_lab) +
    theme(
      panel.grid.major = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  # Add percentage labels if requested
  if (show_percentages) {
    p <- p + geom_text(aes(label = sprintf("%.1f%%", Percent)),
                       position = position_stack(vjust = 0.5),
                       size = text_size, color = "black")
  }
  
  return(p)
}
