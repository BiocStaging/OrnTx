utils::globalVariables(c("log2FoldChange", "neg_log10_p",
                          "color_group", "gene_id"))

#' Differential Expression Analysis
#'
#' @param se A SummarizedExperiment object
#' @param design Formula for experimental design
#' @param contrast Character vector of length 3
#' @param alpha Significance threshold
#' @param lfc_threshold Minimum log2 fold change
#'
#' @return A data.frame of DE results
#' @export
#'
#' @examples
#' data(rosa_example)
#' results <- runDE(rosa_example,
#'   design   = ~ tissue_simple,
#'   contrast = c("tissue_simple", "petal", "abscission zone"))
runDE <- function(se,
                  design        = ~ tissue_clean,
                  contrast      = NULL,
                  alpha         = 0.05,
                  lfc_threshold = 1) {

  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("Install DESeq2: BiocManager::install('DESeq2')")

  message("[OrnTx] Running DESeq2 | design: ", deparse(design))

  # Build DESeq2 object
  dds  <- DESeq2::DESeqDataSet(se, design = design)
  keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
  dds  <- dds[keep, ]
  message("[OrnTx] Genes after filtering: ", nrow(dds))

  # Run DESeq2
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  # Get results
  if (!is.null(contrast)) {
    res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  } else {
    res <- DESeq2::results(dds, alpha = alpha)
  }

  # Convert to data.frame
  res_df             <- as.data.frame(res)
  res_df$gene_id     <- rownames(res_df)
  res_df$significant <- !is.na(res_df$padj) &
    res_df$padj < alpha &
    abs(res_df$log2FoldChange) >= lfc_threshold
  res_df$direction   <- ifelse(!res_df$significant, "NS",
                               ifelse(res_df$log2FoldChange > 0, "Up", "Down"))

  n_up   <- sum(res_df$direction == "Up",   na.rm = TRUE)
  n_down <- sum(res_df$direction == "Down", na.rm = TRUE)
  message("[OrnTx] Up: ", n_up, " | Down: ", n_down,
          " | Total sig: ", n_up + n_down)

  res_df[order(res_df$padj, na.last = TRUE), ]
}


#' Volcano plot of DE results
#'
#' @param de_results Data.frame from runDE()
#' @param alpha Significance threshold
#' @param lfc_threshold LFC threshold
#' @param top_n Number of genes to label
#' @param title Plot title
#'
#' @return A ggplot2 volcano plot
#' @export
#'
#' @examples
#' data(rosa_example)
#' results <- runDE(rosa_example, design = ~ tissue_simple,
#'   contrast = c('tissue_simple', 'petal', 'abscission zone'))
#' plotVolcano(results, title = 'Petal vs AZ')
#' \donttest{
#' data(rosa_example)
#' results <- runDE(rosa_example,
#'   design   = ~ tissue_simple,
#'   contrast = c("tissue_simple", "petal", "abscission zone"))
#' plotVolcano(results, title = "Petal vs Abscission Zone")
#' }
plotVolcano <- function(de_results,
                        alpha         = 0.05,
                        lfc_threshold = 1,
                        top_n         = 20,
                        title         = "Volcano Plot") {
  df              <- de_results
  df$neg_log10_p  <- -log10(pmax(df$padj, 1e-300))
  df$color_group  <- "NS"
  df$color_group[!is.na(df$padj) & df$padj < alpha &
                   df$log2FoldChange >= lfc_threshold]  <- "Up"
  df$color_group[!is.na(df$padj) & df$padj < alpha &
                   df$log2FoldChange <= -lfc_threshold] <- "Down"

  col_vals  <- c("Up" = "#C94040", "Down" = "#3A75C4", "NS" = "#CCCCCC")
  sig_df    <- df[df$color_group != "NS", ]
  top_genes <- utils::head(sig_df[order(-sig_df$neg_log10_p), ], top_n)

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = log2FoldChange,
                                    y = neg_log10_p,
                                    color = color_group)) +
    ggplot2::geom_point(size = 1.5, alpha = 0.6) +
    ggplot2::scale_color_manual(values = col_vals) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", color = "grey40",
                        linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(alpha),
                        linetype = "dashed", color = "grey40",
                        linewidth = 0.4) +
    ggplot2::labs(
      title = title,
      x     = expression(log[2]~"Fold Change"),
      y     = expression(-log[10]~"adjusted p-value"),
      color = NULL) +
    ornTheme()

  if (nrow(top_genes) > 0 && "gene_id" %in% names(top_genes)) {
    p <- p + ggplot2::geom_text(
      data  = top_genes,
      ggplot2::aes(label = gene_id),
      size  = 2.5, color = "grey20", nudge_y = 0.5)
  }
  p
}
