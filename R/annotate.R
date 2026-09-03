utils::globalVariables(c("rosa_annotation", "ann", "gene_ann"))

#' Annotate gene IDs with Rosa chinensis locus information
#'
#' @description
#' Maps generic gene IDs (gene0, gene1, etc.) from the RchiOBHm-V2
#' reference genome to official Rosa chinensis locus names and
#' gene biotype information.
#'
#' @param gene_ids Character vector of gene IDs to annotate
#'   (e.g., c("gene24707", "gene47555")).
#' @param species Character. Currently supports "Rosa chinensis".
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{gene_id} — input gene ID
#'     \item \code{locus_name} — official RchiOBHm-V2 locus name
#'     \item \code{gene_biotype} — gene type (protein_coding, lncRNA, etc.)
#'   }
#' @export
#'
#' @examples
#' ann <- annotateGenes(c("gene10715", "gene21520", "gene37167"))
#' print(ann)
annotateGenes <- function(gene_ids,
                          species = "Rosa chinensis") {

  if (species != "Rosa chinensis") {
    message("[OrnTx] Note: built-in annotation only available ",
            "for Rosa chinensis. Returning gene_ids only.")
    return(data.frame(gene_id      = gene_ids,
                      locus_name   = NA_character_,
                      gene_biotype = NA_character_))
  }

  # Load built-in annotation
  ann_env <- new.env()
  utils::data("rosa_annotation", package = "OrnTx",
       envir = ann_env)
  ann <- ann_env$rosa_annotation

  # Merge
  merged <- merge(
    data.frame(gene_id = gene_ids, stringsAsFactors = FALSE),
    ann,
    by     = "gene_id",
    all.x  = TRUE,
    sort   = FALSE
  )

  # Preserve input order
  merged <- merged[match(gene_ids, merged$gene_id), ]
  rownames(merged) <- NULL

  n_annotated <- sum(!is.na(merged$locus_name))
  message("[OrnTx] Annotated ", n_annotated, " of ",
          length(gene_ids), " genes")
  merged
}


#' Annotate DE results with Rosa chinensis gene information
#'
#' @description
#' Convenience wrapper that adds locus names to a complete DE
#' results data.frame from \code{runDE()}.
#'
#' @param de_results A data.frame from \code{runDE()}.
#' @param species Character. Currently "Rosa chinensis".
#'
#' @return The input data.frame with additional columns:
#'   \code{locus_name} and \code{gene_biotype}.
#' @export
#'
#' @examples
#' \donttest{
#' data(rosa_example)
#' results <- runDE(rosa_example,
#'   design   = ~ tissue_simple,
#'   contrast = c("tissue_simple", "petal", "abscission zone"))
#' results_ann <- annotateDEResults(results)
#' head(results_ann[results_ann$significant == TRUE, ])
#' }
annotateDEResults <- function(de_results,
                              species = "Rosa chinensis") {

  if (!"gene_id" %in% names(de_results))
    stop("de_results must have a 'gene_id' column. ",
         "Make sure you used runDE() from OrnTx.")

  gene_ann <- annotateGenes(de_results$gene_id, species = species)

  result <- cbind(de_results,
                  locus_name   = gene_ann$locus_name,
                  gene_biotype = gene_ann$gene_biotype)
  result
}


#' Get Rosa chinensis annotation database
#'
#' @description
#' Returns the complete built-in annotation table for
#' Rosa chinensis RchiOBHm-V2 genome containing 51,302 genes.
#'
#' @return A data.frame with columns gene_id, locus_name, gene_biotype.
#' @export
#'
#' @examples
#' ann <- getRosaAnnotation()
#' head(ann)
#' nrow(ann)
getRosaAnnotation <- function() {
  ann_env <- new.env()
  utils::data("rosa_annotation", package = "OrnTx", envir = ann_env)
  ann_env$rosa_annotation
}

#' Annotate gene IDs using NCBI functional annotation
#'
#' @description
#' Adds functional product descriptions to DE results when gene IDs
#' are in NCBI LOC format (e.g., LOC112190511). This complements
#' annotateGenes() which works with RchiOBHm locus names.
#'
#' @param de_results A data.frame with a gene_id column containing
#'   NCBI LOC gene identifiers.
#'
#' @return The input data.frame with product column added.
#' @export
#' @name annotateWithNCBI
#'
#' @examples
#' \donttest{
#' data(rosa_ncbi_annotation)
#' de <- data.frame(gene_id = rose_ncbi_annotation$ncbi_gene_id[1:5],
#'                  log2FoldChange = rnorm(5))
#' annotateWithNCBI(de)
#' }
annotateWithNCBI <- function(de_results) {
  ann_env <- new.env()
  utils::data("rose_ncbi_annotation", package = "OrnTx",
              envir = ann_env)
  ann <- ann_env$rose_ncbi_annotation

  merged <- merge(de_results,
                  ann,
                  by.x = "gene_id",
                  by.y = "ncbi_gene_id",
                  all.x = TRUE)

  n_ann <- sum(!is.na(merged$product))
  message("[OrnTx] Annotated ", n_ann, " of ",
          nrow(de_results), " genes with NCBI product names")
  merged
}
