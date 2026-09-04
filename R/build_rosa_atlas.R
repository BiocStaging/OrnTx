utils::globalVariables(c("tissue_clean"))

#' Build Rosa chinensis Expression Atlas
#'
#' @param counts_dir Path to directory containing ReadsPerGene.out.tab files
#' @param metadata_file Path to SraRunTable.txt
#' @param strandedness One of "reverse", "forward", "unstranded"
#' @param min_samples Minimum samples a gene must be expressed in
#'
#' @return A SummarizedExperiment object
#' @export
#'
#' @examples
#' ann <- getRosaAnnotation()
#' nrow(ann)
#' \dontrun{
#' se <- buildRosaAtlas(
#'   counts_dir    = "path/to/counts",
#'   metadata_file = "path/to/SraRunTable.txt"
#' )
#' }
buildRosaAtlas <- function(counts_dir,
                           metadata_file,
                           strandedness = "reverse",
                           min_samples  = 10) {

  message("[OrnTx] Building Rosa chinensis expression atlas...")

  # Step 1: Read metadata
  message("[OrnTx] Reading sample metadata...")
  meta <- utils::read.csv(metadata_file,
                          stringsAsFactors = FALSE,
                          check.names      = FALSE)

  # Keep key columns
  key_cols <- c("Run", "tissue", "Cultivar", "dev_stage",
                "treatment", "BioProject", "geo_loc_name",
                "plant_structure", "Instrument", "LibraryLayout")
  key_cols    <- key_cols[key_cols %in% names(meta)]
  meta_sub    <- meta[, key_cols, drop = FALSE]
  rownames(meta_sub) <- meta_sub$Run

  # Step 2: Find count files
  message("[OrnTx] Scanning count files...")
  count_files <- list.files(counts_dir,
                            pattern    = "_ReadsPerGene.out.tab",
                            full.names = TRUE)
  srr_ids <- gsub("_ReadsPerGene.out.tab", "",
                  basename(count_files))
  names(count_files) <- srr_ids

  # Match with metadata
  common_srr  <- intersect(srr_ids, rownames(meta_sub))
  message("[OrnTx] Samples with counts + metadata: ",
          length(common_srr))

  count_files <- count_files[common_srr]
  meta_final  <- meta_sub[common_srr, , drop = FALSE]

  # Step 3: Choose strand column
  col_idx <- switch(strandedness,
                    reverse    = 4L,
                    forward    = 3L,
                    unstranded = 2L
  )

  # Step 4: Read all count files
  message("[OrnTx] Reading ", length(count_files),
          " count files...")
  count_list <- lapply(seq_along(count_files), function(i) {
    fp  <- count_files[i]
    dat <- utils::read.table(fp, header = FALSE,
                             sep = "\t", skip = 4,
                             stringsAsFactors = FALSE)
    v        <- dat[, col_idx]
    names(v) <- dat[, 1]
    v
  })
  names(count_list) <- names(count_files)

  # Step 5: Build count matrix
  message("[OrnTx] Building count matrix...")
  common_genes <- Reduce(intersect, lapply(count_list, names))
  message("[OrnTx] Common genes: ", length(common_genes))

  count_mat <- do.call(cbind,
                       lapply(count_list, function(x) x[common_genes]))
  storage.mode(count_mat) <- "integer"

  # Step 6: Filter lowly expressed genes
  keep      <- rowSums(count_mat >= 5) >= min_samples
  count_mat <- count_mat[keep, , drop = FALSE]
  message("[OrnTx] Genes after filtering: ", nrow(count_mat))

  # Step 7: Clean tissue metadata
  meta_final$tissue_clean <- tolower(trimws(meta_final$tissue))
  meta_final$tissue_clean <- gsub("petals$", "petal",
                                  meta_final$tissue_clean)
  meta_final$tissue_clean <- ifelse(
    meta_final$tissue_clean == "",
    "unknown",
    meta_final$tissue_clean
  )

  # Step 8: Build SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = count_mat),
    colData = S4Vectors::DataFrame(meta_final)
  )

  message("[OrnTx] Atlas complete: ",
          nrow(se), " genes x ",
          ncol(se), " samples")
  message("[OrnTx] Tissues found: ",
          paste(unique(meta_final$tissue_clean)[1:5],
                collapse = ", "), "...")
  se
}
