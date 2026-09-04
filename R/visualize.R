utils::globalVariables(c("PC1", "PC2", "group"))

#' PCA plot of expression data
#'
#' @param se A SummarizedExperiment object
#' @param assay_name Which assay to use
#' @param color_by Column in colData for point color
#' @param shape_by Column in colData for point shape
#' @param label_points Show sample labels
#'
#' @return A ggplot2 PCA plot
#' @export
#'
#' @examples
#' data(rosa_example)
#' se_norm <- normalizeCounts(rosa_example, method = "CPM")
plotPCA <- function(se,
                    assay_name   = "counts",
                    color_by     = "stage",
                    shape_by     = NULL,
                    label_points = TRUE) {

  expr_mat <- SummarizedExperiment::assay(se, assay_name)
  expr_log <- log2(expr_mat + 1)

  gene_var  <- apply(expr_log, 1, stats::var)
  top_genes <- names(sort(gene_var,
                          decreasing = TRUE))[1:min(500, nrow(expr_log))]
  expr_sub  <- t(expr_log[top_genes, ])

  pca_obj <- stats::prcomp(expr_sub, center = TRUE, scale. = FALSE)
  var_exp <- round(pca_obj$sdev^2 /
                   sum(pca_obj$sdev^2) * 100, 1)

  pca_df <- data.frame(
    PC1    = pca_obj$x[, 1],
    PC2    = pca_obj$x[, 2],
    sample = rownames(pca_obj$x)
  )

  col_df <- as.data.frame(SummarizedExperiment::colData(se))
  common <- intersect(rownames(pca_df), rownames(col_df))
  if (length(common) > 0) {
    pca_df <- cbind(pca_df[common, ],
                    col_df[common, , drop = FALSE])
  }

  if (color_by %in% names(pca_df)) {
    pca_df$group <- pca_df[[color_by]]
  } else {
    pca_df$group <- "Sample"
  }

  rose_colors <- c(
    "#C94040", "#5B8DB8", "#7EB8A4",
    "#E8A598", "#9B7EC8", "#C4956A",
    "#4E9A5A", "#D4B483"
  )

  shape_aes <- if (!is.null(shape_by) &&
                    shape_by %in% names(pca_df)) {
    pca_df[[shape_by]]
  } else {
    NULL
  }

  p <- ggplot2::ggplot(pca_df,
        ggplot2::aes(
          x     = PC1,
          y     = PC2,
          color = group,
          shape = shape_aes,
          label = sample)) +
    ggplot2::geom_point(size = 4, alpha = 0.85) +
    ggplot2::scale_color_manual(values = rose_colors) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(
                       hjust = 0.5,
                       face  = "bold",
                       color = "#2E4057"),
      legend.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = "PCA - Rosa chinensis Expression",
      x     = paste0("PC1 (", var_exp[1], "% variance)"),
      y     = paste0("PC2 (", var_exp[2], "% variance)"),
      color = color_by,
      shape = shape_by
    )

  if (label_points) {
    p <- p + ggplot2::geom_text(
      nudge_y = 0.5,
      size    = 3,
      color   = "grey30"
    )
  }
  p
}

#' OrnTx ggplot2 theme
#'
#' @param base_size Base font size
#' @return A ggplot2 theme
#' @examples
#' p <- ggplot2::ggplot(data.frame(x=1:5, y=1:5),
#'        ggplot2::aes(x, y)) + ornTheme()
#' print(p)
#' @export
ornTheme <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
  ggplot2::theme(
    plot.title   = ggplot2::element_text(
                     hjust = 0.5,
                     face  = "bold",
                     color = "#2E4057",
                     size  = base_size + 2),
    panel.background = ggplot2::element_rect(
                         fill  = "#FAFAF8",
                         color = NA),
    panel.grid.major = ggplot2::element_line(
                         color     = "#E8E8E0",
                         linewidth = 0.3),
    axis.line    = ggplot2::element_line(
                     color     = "#2E4057",
                     linewidth = 0.5),
    legend.background = ggplot2::element_rect(
                          fill  = "white",
                          color = "#CCCCCC")
  )
}

