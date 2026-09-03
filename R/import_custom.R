#' Import RNA-seq data from any plant species into OrnTx
#'
#' @description
#' A universal import function that accepts STAR ReadsPerGene.out.tab
#' files from any plant species. Unlike \code{buildRosaAtlas()} which
#' is Rosa chinensis specific, this function works for any plant species
#' including ornamental flowers, fruit crops, vegetables, and model plants.
#'
#' @param counts_dir Character. Path to directory containing
#'   ReadsPerGene.out.tab files from STAR alignment.
#' @param metadata_file Character. Path to a CSV or TSV file with sample
#'   metadata. Must have a column matching sample IDs. See Details.
#' @param species Character. Scientific name of the plant species,
#'   e.g., \code{"Chrysanthemum morifolium"}, \code{"Tagetes erecta"}.
#' @param sample_col Character. Column name in metadata containing
#'   sample IDs that match count file names (default \code{"Run"}).
#' @param strandedness Character. Library strandedness:
#'   \code{"reverse"} (most dUTP libraries, default),
#'   \code{"forward"}, or \code{"unstranded"}.
#' @param min_count Integer. Minimum count threshold for filtering
#'   (default 5).
#' @param min_samples Integer. Minimum number of samples a gene must
#'   be expressed in (default 3).
#' @param annotation_file Character. Optional path to GTF/GFF file
#'   for gene name annotation. If provided, gene IDs are mapped to
#'   gene names automatically.
#' @param sep Character. Metadata file separator:
#'   \code{","} for CSV (default) or \code{"\t"} for TSV.
#'
#' @return A \code{SummarizedExperiment} object with:
#'   \itemize{
#'     \item \code{assays}: raw counts matrix
#'     \item \code{colData}: sample metadata
#'     \item \code{rowData}: gene annotation (if GTF provided)
#'     \item \code{metadata}: species name and import parameters
#'   }
#'
#' @details
#' The metadata file should contain one row per sample with at minimum
#' a column of sample IDs matching the count file names. For example,
#' if count files are named \code{Sample1_ReadsPerGene.out.tab} and
#' \code{Sample2_ReadsPerGene.out.tab}, the metadata must have a column
#' containing \code{Sample1} and \code{Sample2}.
#'
#' Recommended metadata columns for ornamental plant studies:
#' \itemize{
#'   \item \code{tissue} - tissue type (petal, leaf, root, etc.)
#'   \item \code{dev_stage} - developmental stage
#'   \item \code{cultivar} - cultivar or accession name
#'   \item \code{treatment} - experimental treatment
#'   \item \code{replicate} - biological replicate number
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Chrysanthemum morifolium
#' se <- importCustomAtlas(
#'   counts_dir    = "path/to/chrysanthemum/counts/",
#'   metadata_file = "chrysanthemum_metadata.csv",
#'   species       = "Chrysanthemum morifolium",
#'   strandedness  = "reverse"
#' )
#'
#' # Example 2: Tagetes erecta (marigold)
#' se <- importCustomAtlas(
#'   counts_dir    = "path/to/marigold/counts/",
#'   metadata_file = "marigold_metadata.csv",
#'   species       = "Tagetes erecta",
#'   annotation_file = "marigold_genome.gtf"
#' )
#'
#' # Example 3: Any plant species
#' se <- importCustomAtlas(
#'   counts_dir    = "path/to/counts/",
#'   metadata_file = "metadata.csv",
#'   species       = "Tulipa gesneriana"
#' )
#'
#' # After import - standard analysis works identically
#' se_norm <- normalizeCounts(se, method = "CPM")
#' plotPCA(se_norm, color_by = "tissue")
#' runDE(se, design = ~ tissue,
#'       contrast = c("tissue", "petal", "leaf"))
#' }
importCustomAtlas <- function(counts_dir,
                              metadata_file,
                              species         = "Unknown species",
                              sample_col      = NULL,
                              strandedness    = c("reverse",
                                                  "forward",
                                                  "unstranded"),
                              min_count       = 5L,
                              min_samples     = 3L,
                              annotation_file = NULL,
                              sep             = ",") {

  strandedness <- match.arg(strandedness)

  message("[OrnTx] Importing data for: ", species)
  message("[OrnTx] Strandedness: ", strandedness)

  # -- Step 1: Read metadata ----------------------------------------------
  message("[OrnTx] Reading metadata from: ",
          basename(metadata_file))

  if (!file.exists(metadata_file))
    stop("Metadata file not found: ", metadata_file)

  meta <- utils::read.table(metadata_file,
                            header           = TRUE,
                            sep              = sep,
                            stringsAsFactors = FALSE,
                            check.names      = FALSE,
                            quote            = "\"",
                            fill             = TRUE)

  message("[OrnTx] Metadata: ", nrow(meta), " rows x ",
          ncol(meta), " columns")

  # Auto-detect sample ID column
  if (is.null(sample_col)) {
    candidate_cols <- c("Run", "SampleID", "Sample", "sample_id",
                        "SRR", "sample", "ID", "id", "name")
    found <- candidate_cols[candidate_cols %in% names(meta)]
    if (length(found) == 0) {
      message("[OrnTx] Could not auto-detect sample ID column.")
      message("[OrnTx] Available columns: ",
              paste(names(meta)[seq_len(min(10, ncol(meta)))], collapse = ", "))
      stop("Please specify 'sample_col' argument with the column ",
           "name containing sample IDs.")
    }
    sample_col <- found[1]
    message("[OrnTx] Auto-detected sample ID column: '",
            sample_col, "'")
  }

  if (!sample_col %in% names(meta))
    stop("Column '", sample_col, "' not found in metadata. ",
         "Available: ", paste(names(meta)[1:10], collapse = ", "))

  rownames(meta) <- make.unique(as.character(meta[[sample_col]]))

  # -- Step 2: Find count files -------------------------------------------
  message("[OrnTx] Scanning: ", counts_dir)

  if (!dir.exists(counts_dir))
    stop("counts_dir not found: ", counts_dir)

  count_files <- list.files(counts_dir,
                            pattern    = "ReadsPerGene.out.tab",
                            full.names = TRUE,
                            recursive  = TRUE)

  if (length(count_files) == 0)
    stop("No ReadsPerGene.out.tab files found in: ", counts_dir,
         "\nMake sure STAR was run with --quantMode GeneCounts")

  # Extract sample IDs from filenames
  sample_ids <- gsub("_ReadsPerGene.out.tab", "",
                     basename(count_files))
  sample_ids <- gsub("\\.ReadsPerGene\\.out\\.tab", "",
                     sample_ids)
  names(count_files) <- sample_ids

  message("[OrnTx] Found ", length(count_files), " count files")

  # -- Step 3: Match samples with metadata -------------------------------
  common <- intersect(sample_ids, rownames(meta))

  if (length(common) == 0) {
    message("[OrnTx] ERROR: No overlap between count files and metadata!")
    message("[OrnTx] Count file IDs (first 5): ",
            paste(sample_ids[seq_len(min(5, length(sample_ids)))], collapse = ", "))
    message("[OrnTx] Metadata IDs (first 5): ",
            paste(rownames(meta)[seq_len(min(5, nrow(meta)))], collapse = ", "))
    stop("Sample IDs in count files don't match metadata. ",
         "Check the 'sample_col' argument.")
  }

  if (length(common) < length(sample_ids)) {
    message("[OrnTx] Note: ", length(sample_ids) - length(common),
            " count files have no matching metadata - skipping those.")
  }

  count_files <- count_files[common]
  meta_final  <- meta[common, , drop = FALSE]
  message("[OrnTx] Matched samples: ", length(common))

  # -- Step 4: Choose strand column --------------------------------------
  col_idx <- switch(strandedness,
                    reverse    = 4L,
                    forward    = 3L,
                    unstranded = 2L
  )

  # -- Step 5: Read count files -------------------------------------------
  message("[OrnTx] Reading ", length(count_files),
          " count files...")

  count_list <- lapply(seq_along(count_files), function(i) {
    fp <- count_files[i]

    # Try reading - catch errors per file
    tryCatch({
      dat <- utils::read.table(fp,
                               header           = FALSE,
                               sep              = "\t",
                               skip             = 4,
                               stringsAsFactors = FALSE)
      if (ncol(dat) < col_idx)
        stop("File has fewer columns than expected for strandedness '",
             strandedness, "': ", basename(fp))
      v        <- dat[, col_idx]
      names(v) <- dat[, 1]
      v
    }, error = function(e) {
      warning("Failed to read: ", basename(fp), " - ", e$message)
      NULL
    })
  })
  names(count_list) <- names(count_files)

  # Remove failed files
  failed <- vapply(count_list, is.null, logical(1))
  if (any(failed)) {
    message("[OrnTx] Warning: ", sum(failed),
            " files failed to read and were excluded.")
    count_list  <- count_list[!failed]
    meta_final  <- meta_final[names(count_list), , drop = FALSE]
  }

  if (length(count_list) == 0)
    stop("No count files could be read successfully.")

  # -- Step 6: Build count matrix -----------------------------------------
  message("[OrnTx] Building count matrix...")
  common_genes <- Reduce(intersect, lapply(count_list, names))
  message("[OrnTx] Common genes across all samples: ",
          length(common_genes))

  if (length(common_genes) == 0)
    stop("No common genes found across samples. ",
         "Check that all files use the same reference genome.")

  count_mat <- do.call(cbind,
                       lapply(count_list, function(x) x[common_genes]))
  colnames(count_mat) <- names(count_list)
  storage.mode(count_mat) <- "integer"

  # -- Step 7: Filter lowly expressed genes ------------------------------
  keep      <- rowSums(count_mat >= min_count) >= min_samples
  count_mat <- count_mat[keep, , drop = FALSE]
  message("[OrnTx] Genes after filtering (>= ", min_count,
          " counts in >= ", min_samples, " samples): ",
          nrow(count_mat))

  # -- Step 8: Gene annotation (optional) --------------------------------
  row_data <- S4Vectors::DataFrame(gene_id = rownames(count_mat))
  rownames(row_data) <- rownames(count_mat)

  if (!is.null(annotation_file)) {
    message("[OrnTx] Parsing gene annotation from: ",
            basename(annotation_file))
    ann <- tryCatch(
      .parseGTFAnnotation(annotation_file,
                          rownames(count_mat)),
      error = function(e) {
        warning("Annotation parsing failed: ", e$message,
                "\nProceeding without annotation.")
        NULL
      }
    )
    if (!is.null(ann)) {
      row_data <- S4Vectors::DataFrame(ann)
      message("[OrnTx] Annotated ", sum(!is.na(ann$gene_name)),
              " genes with names")
    }
  }

  # -- Step 9: Build SummarizedExperiment --------------------------------
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = count_mat),
    colData = S4Vectors::DataFrame(meta_final),
    rowData = row_data,
    metadata = list(
      species      = species,
      strandedness = strandedness,
      n_samples    = ncol(count_mat),
      n_genes      = nrow(count_mat),
      import_date  = Sys.Date(),
      ornAtlas_version = utils::packageVersion("OrnTx")
    )
  )

  message("[OrnTx] OK Import complete!")
  message("[OrnTx] Species:  ", species)
  message("[OrnTx] Genes:    ", nrow(se))
  message("[OrnTx] Samples:  ", ncol(se))
  message("[OrnTx] Ready for: normalizeCounts(), plotPCA(), runDE()")
  se
}


# -- Internal: Parse GTF annotation ----------------------------------------

#' @keywords internal
.parseGTFAnnotation <- function(gtf_file, gene_ids) {

  if (!file.exists(gtf_file))
    stop("GTF file not found: ", gtf_file)

  # Read GTF - only gene lines for speed
  message("[OrnTx] Reading GTF (this may take a moment)...")
  gtf <- utils::read.table(gtf_file,
                           header           = FALSE,
                           sep              = "\t",
                           comment.char     = "#",
                           stringsAsFactors = FALSE,
                           quote            = "")

  # Keep only 'gene' feature lines
  gene_lines <- gtf[gtf[, 3] == "gene", , drop = FALSE]
  if (nrow(gene_lines) == 0)
    gene_lines <- gtf  # fallback if no 'gene' lines

  message("[OrnTx] Parsing ", nrow(gene_lines), " gene records...")

  # Extract gene_id and gene_name from attributes column
  attrs      <- gene_lines[, 9]
  gene_id    <- .extractGTFAttr(attrs, "gene_id")
  gene_name  <- .extractGTFAttr(attrs, "gene_name")
  gene_biotype <- .extractGTFAttr(attrs, "gene_biotype")

  ann_df <- data.frame(
    gene_id      = gene_id,
    gene_name    = gene_name,
    gene_biotype = gene_biotype,
    stringsAsFactors = FALSE
  )
  ann_df <- ann_df[!duplicated(ann_df$gene_id), ]
  rownames(ann_df) <- ann_df$gene_id

  # Match to our genes
  matched <- ann_df[intersect(gene_ids, rownames(ann_df)), ,
                    drop = FALSE]

  # Add missing genes with NA
  missing      <- setdiff(gene_ids, rownames(matched))
  if (length(missing) > 0) {
    missing_df <- data.frame(
      gene_id      = missing,
      gene_name    = NA_character_,
      gene_biotype = NA_character_,
      row.names    = missing,
      stringsAsFactors = FALSE
    )
    matched <- rbind(matched, missing_df)
  }
  matched[gene_ids, , drop = FALSE]
}


#' @keywords internal
.extractGTFAttr <- function(attrs, attr_name) {
  pattern <- paste0(attr_name, ' "([^"]+)"')
  matches <- regmatches(attrs, regexpr(pattern, attrs))
  result  <- rep(NA_character_, length(attrs))
  has_match <- nchar(matches) > 0
  if (any(has_match)) {
    extracted <- gsub(paste0(attr_name, ' "([^"]+)"'), "\\1",
                      matches[has_match])
    result[has_match] <- extracted
  }
  result
}
