utils::globalVariables(c("gene", "stage", "expression"))

#' Build an expression atlas from count data
#'
#' @param se A SummarizedExperiment object from importRosaCounts()
#' @param group_var Column name in colData for grouping samples
#' @param assay_name Which assay to use (default "counts")
#'
#' @return A SummarizedExperiment with atlas slot
#' @export
#'
#' @examples
#' library(OrnTx)
#' data(rosa_example)
#' se_norm <- normalizeCounts(rosa_example, method = "CPM")
buildAtlas <- function(se,
                       group_var  = "stage",
                       assay_name = "counts") {

  # Get expression matrix
  expr_mat <- SummarizedExperiment::assay(se, assay_name)

  # Get sample groups
  col_df <- as.data.frame(SummarizedExperiment::colData(se))

  if (!group_var %in% names(col_df))
    stop("'", group_var, "' not found in sample metadata.")

  groups <- col_df[[group_var]]
  message("Building atlas for groups: ",
          paste(unique(groups), collapse = ", "))

  # Average expression per group
  group_levels <- unique(groups)
  atlas_mat <- vapply(group_levels, function(g) {
    idx <- which(groups == g)
    if (length(idx) == 1) expr_mat[, idx]
    else rowMeans(expr_mat[, idx, drop = FALSE])
  }, FUN.VALUE = numeric(nrow(expr_mat)))
  rownames(atlas_mat) <- rownames(expr_mat)

  # Filter lowly expressed genes
  keep <- rowSums(atlas_mat >= 10) >= 1
  atlas_mat <- atlas_mat[keep, , drop = FALSE]
  message("Genes retained after filtering: ", nrow(atlas_mat))

  # Build result object
  atlas_se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(expression = atlas_mat),
    colData = S4Vectors::DataFrame(
      stage = group_levels
    )
  )
  atlas_se
}


#' Plot expression atlas heatmap
#'
#' @param atlas A SummarizedExperiment from buildAtlas()
#' @param top_n Number of top variable genes to show
#' @param title Plot title
#'
#' @return A ggplot2 heatmap object
#' @export
#' @examples
#' data(rosa_example)
#' se_norm <- normalizeCounts(rosa_example, method = 'CPM')
#' atlas <- buildAtlas(se_norm, group_var = 'tissue_simple',
#'                     assay_name = 'CPM')
#' plotAtlas(atlas, top_n = 30, title = 'Rosa chinensis Atlas')
#' \donttest{
#' data(rosa_example)
#' se_norm <- normalizeCounts(rosa_example, method = "CPM")
#' atlas <- buildAtlas(se_norm, group_var = "tissue_simple",
#'                     assay_name = "CPM")
#' plotAtlas(atlas, top_n = 30, title = "Rosa chinensis Atlas")
#' }

plotAtlas <- function(atlas,
                      top_n = 50,
                      title = "Expression Atlas") {
  expr_mat <- SummarizedExperiment::assay(atlas, "expression")

  # Select top variable genes
  gene_var  <- apply(expr_mat, 1, stats::var)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(top_n, nrow(expr_mat))]
  plot_mat  <- expr_mat[top_genes, , drop = FALSE]

  # Scale per gene (z-score)
  plot_mat <- t(scale(t(plot_mat)))

  # Convert to long format for ggplot
  df <- as.data.frame(plot_mat)
  df$gene  <- rownames(df)
  df_long  <- tidyr::pivot_longer(df,
                                  cols      = -gene,
                                  names_to  = "stage",
                                  values_to = "expression")

  # Rose-inspired color palette
  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = stage, y = gene, fill = expression)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(
      low      = "#2E4057",
      mid      = "#FAFAF8",
      high     = "#C94040",
      midpoint = 0,
      name     = "Z-score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.title   = ggplot2::element_text(
        hjust = 0.5, face = "bold",
        color = "#2E4057"),
      panel.grid   = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      x     = "Developmental Stage",
      y     = paste0("Top ", top_n, " Variable Genes")
    )
}
