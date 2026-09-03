
#' Import STAR count files into OrnTx
#'
#' @param file_paths Named character vector of paths to
#'   ReadsPerGene.out.tab files.
#' @param col_data A data.frame with sample metadata.
#' @param strandedness One of "reverse", "forward", "unstranded".
#'
#' @return A SummarizedExperiment object.
#' @export
#'
#' @examples
#' \dontrun{
#' files <- c(Sample1 = "path/to/sample1/ReadsPerGene.out.tab")
#' importRosaCounts(files, col_data = metadata)
#' }
importRosaCounts <- function(file_paths,
                              col_data,
                              strandedness = "reverse") {
  col_idx <- switch(strandedness,
    reverse    = 4L,
    forward    = 3L,
    unstranded = 2L
  )
  count_list <- lapply(seq_along(file_paths), function(i) {
    dat <- utils::read.table(file_paths[i], header = FALSE,
                              sep = "	", skip = 4,
                              stringsAsFactors = FALSE)
    v <- dat[, col_idx]
    names(v) <- dat[, 1]
    v
  })
  names(count_list) <- names(file_paths)
  common_genes <- Reduce(intersect, lapply(count_list, names))
  count_mat <- do.call(cbind,
    lapply(count_list, function(x) x[common_genes]))
  storage.mode(count_mat) <- "integer"
  SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = count_mat),
    colData = S4Vectors::DataFrame(col_data)
  )
}


#' Normalize count data
#'
#' @param se A SummarizedExperiment object.
#' @param method One of "CPM", "TMM", or "DESeq2_vst".
#' @param design Optional formula for DESeq2 methods.
#'
#' @return SummarizedExperiment with normalized assay added.
#' @export
#'
#' @examples
#' data(rosa_example)
#' se_norm <- normalizeCounts(rosa_example, method = "CPM")
normalizeCounts <- function(se,
                             method = c("CPM", "TMM", "DESeq2_vst"),
                             design = ~ 1) {
  method     <- match.arg(method)
  counts_mat <- SummarizedExperiment::assay(se, "counts")
  norm_mat <- switch(method,
    "CPM" = {
      t(t(counts_mat) / (colSums(counts_mat) / 1e6))
    },
    "TMM" = {
      if (!requireNamespace("edgeR", quietly = TRUE))
        stop("Install edgeR: BiocManager::install('edgeR')")
      dge <- edgeR::DGEList(counts = counts_mat)
      dge <- edgeR::calcNormFactors(dge, method = "TMM")
      edgeR::cpm(dge, log = TRUE, prior.count = 1)
    },
    "DESeq2_vst" = {
      if (!requireNamespace("DESeq2", quietly = TRUE))
        stop("Install DESeq2: BiocManager::install('DESeq2')")
      dds <- DESeq2::DESeqDataSet(se, design = design)
      dds <- DESeq2::estimateSizeFactors(dds)
      as.matrix(SummarizedExperiment::assay(
        DESeq2::vst(dds, blind = TRUE)))
    }
  )
  SummarizedExperiment::assay(se, method) <- norm_mat
  message("[OrnTx] Normalization complete: assay '",
          method, "' added.")
  se
}


#' Import count data from any RNA-seq tool
#'
#' @description
#' Accepts count matrices from STAR, featureCounts, HTSeq,
#' Salmon, or any pre-loaded matrix. Unlike importCustomAtlas()
#' which reads STAR output files, this function works with
#' already-loaded matrices or files from any quantification tool.
#'
#' @param counts File path or pre-loaded matrix/data.frame with
#'   genes as rows and samples as columns.
#' @param metadata Data.frame with sample metadata. Row names
#'   or first column must match column names of counts.
#' @param species Scientific name of the plant species.
#' @param file_format One of "auto", "STAR", "featureCounts",
#'   "HTSeq", "Salmon", "matrix". Default "auto" detects format.
#' @param sample_col Column in metadata containing sample IDs.
#' @param min_count Minimum count for filtering (default 5).
#' @param min_samples Minimum samples expressing gene (default 3).
#'
#' @return A SummarizedExperiment object.
#' @export
#'
#' @examples
#' counts <- matrix(rpois(500, 10), nrow=50, ncol=10)
#' rownames(counts) <- paste0("gene", 1:50)
#' colnames(counts) <- paste0("sample", 1:10)
#' meta <- data.frame(Sample = colnames(counts),
#'                    tissue = rep(c("petal","leaf"), 5))
#' se <- importCountMatrix(counts, meta,
#'   species = "Rosa chinensis", sample_col = "Sample")
importCountMatrix <- function(counts,
                               metadata,
                               species     = "Unknown species",
                               file_format = "auto",
                               sample_col  = NULL,
                               min_count   = 5L,
                               min_samples = 3L) {

  message("[OrnTx] Importing count data for: ", species)

  # Read counts from file or use matrix directly
  if (is.character(counts) && length(counts) == 1) {
    if (!file.exists(counts))
      stop("Count file not found: ", counts)

    if (file_format == "auto") {
      first_line <- readLines(counts, n = 1)
      file_format <-
        if (grepl("Geneid.*Chr.*Start", first_line)) "featureCounts"
        else "matrix"
      message("[OrnTx] Auto-detected format: ", file_format)
    }

    count_mat <- switch(file_format,
      "featureCounts" = {
        dat <- utils::read.table(counts, header = TRUE,
                                  sep = "	", skip = 1,
                                  stringsAsFactors = FALSE,
                                  check.names = FALSE)
        mat <- as.matrix(dat[, 7:ncol(dat)])
        rownames(mat) <- dat[, 1]
        colnames(mat) <- gsub(".*/|[.]bam$|[.]BAM$", "",
                               colnames(mat))
        mat
      },
      "HTSeq" = {
        dat <- utils::read.table(counts, header = FALSE,
                                  sep = "	",
                                  stringsAsFactors = FALSE)
        dat <- dat[!grepl("^__", dat[, 1]), ]
        mat <- matrix(dat[, 2], ncol = 1)
        rownames(mat) <- dat[, 1]
        colnames(mat) <- basename(counts)
        mat
      },
      "STAR" = {
        dat <- utils::read.table(counts, header = FALSE,
                                  sep = "	", skip = 4,
                                  stringsAsFactors = FALSE)
        mat <- matrix(dat[, 4], ncol = 1)
        rownames(mat) <- dat[, 1]
        colnames(mat) <- basename(counts)
        mat
      },
      "Salmon" = {
        dat <- utils::read.table(counts, header = TRUE,
                                  sep = "	",
                                  stringsAsFactors = FALSE)
        col <- if ("NumReads" %in% names(dat)) "NumReads"
               else names(dat)[2]
        mat <- matrix(as.integer(round(dat[[col]])), ncol = 1)
        rownames(mat) <- dat[, 1]
        colnames(mat) <- basename(counts)
        mat
      },
      {
        dat <- utils::read.table(counts, header = TRUE,
                                  sep = "	", row.names = 1,
                                  stringsAsFactors = FALSE,
                                  check.names = FALSE)
        as.matrix(dat)
      }
    )
  } else if (is.matrix(counts) || is.data.frame(counts)) {
    count_mat <- as.matrix(counts)
    message("[OrnTx] Using pre-loaded count matrix: ",
            nrow(count_mat), " genes x ",
            ncol(count_mat), " samples")
  } else {
    stop("counts must be a file path, matrix, or data.frame.")
  }

  storage.mode(count_mat) <- "integer"

  # Read metadata
  if (is.character(metadata) && length(metadata) == 1) {
    if (!file.exists(metadata))
      stop("Metadata file not found: ", metadata)
    sep <- if (grepl("[.]tsv$|[.]tab$", metadata)) "	" else ","
    metadata <- utils::read.table(metadata, header = TRUE,
                                   sep = sep,
                                   stringsAsFactors = FALSE,
                                   check.names = FALSE)
  }
  meta <- as.data.frame(metadata)

  # Auto-detect sample column
  if (is.null(sample_col)) {
    candidates <- c("Run","SampleID","Sample","sample","ID","id")
    found <- candidates[candidates %in% names(meta)]
    if (length(found) > 0) {
      sample_col <- found[1]
      message("[OrnTx] Sample ID column: ", sample_col)
    }
  }
  if (!is.null(sample_col) && sample_col %in% names(meta))
    rownames(meta) <- meta[[sample_col]]

  # Match samples
  common <- intersect(colnames(count_mat), rownames(meta))
  if (length(common) == 0)
    stop("No matching samples between counts and metadata.")
  count_mat <- count_mat[, common, drop = FALSE]
  meta      <- meta[common, , drop = FALSE]
  message("[OrnTx] Matched samples: ", length(common))

  # Filter
  keep      <- rowSums(count_mat >= min_count) >= min_samples
  count_mat <- count_mat[keep, , drop = FALSE]
  message("[OrnTx] Genes after filtering: ", nrow(count_mat))

  # Build SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays   = list(counts = count_mat),
    colData  = S4Vectors::DataFrame(meta),
    metadata = list(
      species     = species,
      file_format = file_format,
      n_genes     = nrow(count_mat),
      n_samples   = ncol(count_mat),
      import_date = Sys.Date(),
      ornAtlas_version = utils::packageVersion("OrnTx")
    )
  )

  message("[OrnTx] Import complete!")
  message("[OrnTx] Species:  ", species)
  message("[OrnTx] Genes:    ", nrow(se))
  message("[OrnTx] Samples:  ", ncol(se))
  message("[OrnTx] Ready for: normalizeCounts(), plotPCA(), runDE()")
  se
}

