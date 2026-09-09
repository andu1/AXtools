#' Bar plot of average UMIs per gene by MSigDB pathway
#'
#' Scores each gene set from MSigDB by its mean UMI count per gene across
#' cells, highlights pathways above a threshold, and returns the plot (and
#' optionally the filtered gene-set list for downstream use).
#'
#' @param seurat_obj A Seurat object.
#' @param clusters Optional cluster identities to subset to; `NULL` uses all
#'   cells.
#' @param species Species name passed to [msigdbr::msigdbr()] (default
#'   `"Mus musculus"`).
#' @param category MSigDB collection letter (default `"H"` for hallmark).
#' @param assay Seurat assay to pull counts from.
#' @param umi_threshold Minimum average UMIs per gene to flag a pathway as
#'   expressed.
#' @param top_n If set, show only the top N pathways by UMI count.
#' @param title Plot title; auto-generated if `NULL`.
#' @param fill_color,highlight_color Bar colours for below- and above-threshold
#'   pathways.
#' @param base_size Base font size for the plot theme.
#' @param return_pathways If `TRUE`, return a list with `plot` and `pathways`
#'   (the filtered gene-set list); otherwise return the plot alone.
#' @return A ggplot object, or a list with `plot` and `pathways` when
#'   `return_pathways = TRUE`.
#' @export
plot_pathway_umis <- function(seurat_obj,
                              clusters = NULL,
                              species = "Mus musculus",
                              category = "H",
                              assay = "RNA",
                              umi_threshold = 0.05,
                              top_n = NULL,
                              title = NULL,
                              fill_color = "steelblue",
                              highlight_color = "darkred",
                              base_size = 14,
                              return_pathways = FALSE) {
  
  # Load required packages
  require(msigdbr)
  require(ggplot2)
  require(Matrix)
  require(Seurat)
  
  # Get MSigDB gene sets
  msigdb_data <- msigdbr(species = species, category = category)
  geneset_list <- split(msigdb_data$gene_symbol, msigdb_data$gs_name)
  geneset_list <- lapply(geneset_list, unique)
  
  # Get expression data (updated to use layer)
  counts <- GetAssayData(seurat_obj, assay = assay, layer = "counts")
  
  # Filter gene sets to genes present in data
  geneset_filtered <- lapply(geneset_list, function(genes) {
    genes[genes %in% rownames(counts)]
  })
  geneset_filtered <- geneset_filtered[sapply(geneset_filtered, length) > 0]
  
  # Subset to clusters of interest if specified
  if (!is.null(clusters)) {
    Idents(seurat_obj) <- "seurat_clusters"
    cells_to_use <- WhichCells(seurat_obj, idents = clusters)
    counts_subset <- counts[, cells_to_use, drop = FALSE]
  } else {
    counts_subset <- counts
  }
  
  # Calculate average UMIs per pathway
  avg_pathway_umis <- sapply(geneset_filtered, function(genes) {
    genes_in_data <- intersect(genes, rownames(counts_subset))
    total_umis <- Matrix::colSums(counts_subset[genes_in_data, , drop = FALSE])
    mean(total_umis) / length(genes_in_data)
  })
  
  # Identify good pathways based on UMI threshold
  good_pathways <- names(avg_pathway_umis)[avg_pathway_umis > umi_threshold]
  
  # Create filtered pathway gene list
  good_pathways_genes <- geneset_filtered[good_pathways]
  good_pathways_genes <- lapply(good_pathways_genes, function(genes) {
    intersect(genes, rownames(seurat_obj))
  })
  good_pathways_genes <- good_pathways_genes[sapply(good_pathways_genes, length) > 0]
  
  # Convert to data frame for plotting (include ALL pathways)
  pathway_df <- data.frame(
    Pathway = names(avg_pathway_umis),
    Avg_UMIs_per_gene = as.numeric(avg_pathway_umis),
    Above_Threshold = names(avg_pathway_umis) %in% names(good_pathways_genes),
    stringsAsFactors = FALSE
  )
  
  # Sort by UMI count
  pathway_df <- pathway_df[order(pathway_df$Avg_UMIs_per_gene, decreasing = TRUE), ]
  
  # Filter to top N if specified
  if (!is.null(top_n)) {
    pathway_df <- head(pathway_df, top_n)
  }
  
  # Clean pathway names (remove HALLMARK_ prefix if present) BEFORE setting factor
  pathway_df$Pathway <- gsub("HALLMARK_", "", pathway_df$Pathway)
  pathway_df$Pathway <- gsub("_", " ", pathway_df$Pathway)
  
  # Set factor levels for ordering (AFTER cleaning names and sorting)
  pathway_df$Pathway <- factor(pathway_df$Pathway, levels = pathway_df$Pathway)
  
  # Create title if not provided
  if (is.null(title)) {
    if (!is.null(clusters)) {
      cluster_text <- paste(clusters, collapse = ", ")
      title <- paste0("Average UMIs per Gene by Pathway (Clusters ", cluster_text, ")")
    } else {
      title <- "Average UMIs per Gene by Pathway (All Cells)"
    }
  }
  
  # Add threshold line text
  subtitle <- paste0("Threshold: ", umi_threshold, " UMIs/gene (",
                     sum(pathway_df$Above_Threshold), " pathways above threshold)")
  
  # Create plot with color based on threshold
  p <- ggplot(pathway_df, aes(x = Pathway, y = Avg_UMIs_per_gene, fill = Above_Threshold)) +
    geom_col() +
    geom_hline(yintercept = umi_threshold, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("FALSE" = fill_color, "TRUE" = highlight_color),
                      labels = c("FALSE" = "Below threshold", "TRUE" = "Above threshold"),
                      name = "") +
    coord_flip() +
    labs(title = title, subtitle = subtitle, x = "Pathway", y = "Avg UMIs per Gene") +
    theme_minimal(base_size = base_size) +
    theme(axis.text.y = element_text(size = base_size - 2),
          legend.position = "bottom")
  
  # Return plot and optionally the pathway list
  if (return_pathways) {
    return(list(plot = p, pathways = good_pathways_genes))
  } else {
    return(p)
  }
}
