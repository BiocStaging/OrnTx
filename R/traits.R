utils::globalVariables(c("floriculture_traits"))

#' @title Get curated floriculture trait gene sets
#' @description Returns curated gene sets anchored to Arabidopsis
#'   orthologs for cross-species trait analysis.
#' @param trait Character. One of "flower_color", "scent",
#'   "senescence", "development", "regulatory", or "all".
#' @param species Character. "Rosa chinensis", "Tagetes erecta",
#'   or "all".
#' @param evidence Character. "Strong", "Moderate", "Weak", or "all".
#' @return A data.frame of curated trait genes.
#' @export
#' @name getFlowerGenes
#' @examples
#' color_genes <- getFlowerGenes(trait = "flower_color")
#' head(color_genes)
getFlowerGenes <- function(trait="all", species="all", evidence="all") {
  db_env <- new.env()
  utils::data("floriculture_traits", package="OrnTx", envir=db_env)
  db <- db_env$floriculture_traits
  if (!identical(trait, "all"))
    db <- db[grepl(paste(trait, collapse="|"), db$trait,
                   ignore.case=TRUE), ]
  if (!identical(species, "all") && !identical(evidence, "all")) {
    if (grepl("rosa|rose", species, ignore.case=TRUE))
      db <- db[db$rose_evidence %in% evidence, ]
    else if (grepl("tagetes|marigold", species, ignore.case=TRUE))
      db <- db[db$marigold_evidence %in% evidence, ]
  }
  message("[OrnTx] Trait genes retrieved: ", nrow(db))
  db
}

#' @title Find orthologs in floriculture trait database
#' @description Maps gene symbols or Arabidopsis IDs to curated
#'   trait gene information.
#' @param gene_ids Character vector of gene IDs to query.
#' @param id_type One of "symbol" or "arabidopsis".
#' @param species Source species (informational).
#' @return A data.frame of matched trait genes.
#' @export
#' @name findOrthologs
#' @examples
#' findOrthologs(c("CHS", "DFR", "ANS"), id_type = "symbol")
findOrthologs <- function(gene_ids, id_type=c("symbol","arabidopsis"),
                           species="any") {
  id_type <- match.arg(id_type)
  db_env  <- new.env()
  utils::data("floriculture_traits", package="OrnTx", envir=db_env)
  db <- db_env$floriculture_traits
  matched <- if (id_type == "symbol") db[db$gene_symbol %in% gene_ids, ]
             else db[db$arabidopsis_id %in% gene_ids, ]
  message("[OrnTx] Queried: ", length(gene_ids),
          " | Matched: ", nrow(matched))
  matched
}

#' @title Annotate DE results with trait information
#' @description Identifies known floriculture trait genes in
#'   differential expression results.
#' @param de_results A data.frame from runDE().
#' @param trait Trait categories to check, or 'all'.
#' @param species Species for evidence filtering, or 'all'.
#' @param gene_col Optional column name containing gene symbols.
#' @param mapping Optional data.frame with gene_id and
#'   arabidopsis_id columns for ortholog-based matching.
#' @export
#' @name annotateWithTraits
#' @examples
#' \donttest{
#' data(rosa_example)
#' res <- runDE(rosa_example, design=~tissue_simple,
#'   contrast=c("tissue_simple","petal","abscission zone"))
#' res_ann <- annotateDEResults(res)
#' annotateWithTraits(res_ann)
#' }
annotateWithTraits <- function(de_results,
                                trait       = "all",
                                species     = "all",
                                gene_col    = NULL,
                                mapping     = NULL) {
  trait_db <- getFlowerGenes(trait=trait, species=species)

  de_results$trait_category <- NA_character_
  de_results$pathway        <- NA_character_
  de_results$arabidopsis_id <- NA_character_
  de_results$is_trait_gene  <- FALSE

  # Strategy 1: User provides a mapping table
  # mapping should have columns: gene_id + gene_symbol or arabidopsis_id
  if (!is.null(mapping)) {
    if ("arabidopsis_id" %in% names(mapping)) {
      for (i in seq_len(nrow(trait_db))) {
        at_id <- trait_db$arabidopsis_id[i]
        hit_map <- mapping$gene_id[mapping$arabidopsis_id == at_id]
        hit <- de_results$gene_id %in% hit_map
        if (any(hit, na.rm=TRUE)) {
          de_results$trait_category[hit] <- trait_db$trait_category[i]
          de_results$pathway[hit]        <- trait_db$pathway[i]
          de_results$arabidopsis_id[hit] <- at_id
          de_results$is_trait_gene[hit]  <- TRUE
        }
      }
    }
  # Strategy 2: User specifies a column with gene symbols
  } else if (!is.null(gene_col) && gene_col %in% names(de_results)) {
    for (i in seq_len(nrow(trait_db))) {
      sym <- trait_db$gene_symbol[i]
      hit <- de_results[[gene_col]] == sym
      hit[is.na(hit)] <- FALSE
      if (any(hit)) {
        de_results$trait_category[hit] <- trait_db$trait_category[i]
        de_results$pathway[hit]        <- trait_db$pathway[i]
        de_results$arabidopsis_id[hit] <- trait_db$arabidopsis_id[i]
        de_results$is_trait_gene[hit]  <- TRUE
      }
    }
  # Strategy 3: Match gene symbols directly in gene_id column
  } else if ("gene_id" %in% names(de_results)) {
    for (i in seq_len(nrow(trait_db))) {
      sym <- trait_db$gene_symbol[i]
      hit <- de_results$gene_id == sym
      hit[is.na(hit)] <- FALSE
      if (any(hit)) {
        de_results$trait_category[hit] <- trait_db$trait_category[i]
        de_results$pathway[hit]        <- trait_db$pathway[i]
        de_results$arabidopsis_id[hit] <- trait_db$arabidopsis_id[i]
        de_results$is_trait_gene[hit]  <- TRUE
      }
    }
  }

  n_trait     <- sum(de_results$is_trait_gene, na.rm=TRUE)
  n_sig_trait <- sum(de_results$is_trait_gene &
                     de_results$significant, na.rm=TRUE)
  message("[OrnTx] Trait genes: ", n_trait)
  message("[OrnTx] Significant trait genes: ", n_sig_trait)
  de_results
}
summarizeTraitExpression <- function(se, trait="flower_color_anthocyanin",
                                      assay_name="CPM", species="all") {
  trait_db <- getFlowerGenes(trait=trait, species=species)
  expr_mat <- SummarizedExperiment::assay(se, assay_name)
  matched  <- rownames(expr_mat)[
    rownames(expr_mat) %in% trait_db$gene_symbol |
    rownames(expr_mat) %in% trait_db$arabidopsis_id]
  if (length(matched) == 0) {
    message("[OrnTx] No trait genes found. Use annotateWithTraits()",
            " on DE results instead.")
    return(invisible(NULL))
  }
  message("[OrnTx] Trait genes found: ", length(matched))
  as.data.frame(expr_mat[matched, , drop=FALSE])
}

#' @title Standardize developmental stage vocabulary
#' @description Harmonizes inconsistent developmental stage names
#'   across studies into a controlled vocabulary. Solves the problem
#'   of "bud", "Bud", "floral bud", "Stage 1" all meaning the same
#'   thing in different publications.
#' @param stage_vector Character vector of stage names to standardize.
#' @return Character vector of standardized stage names.
#' @export
#' @name standardizeStage
#' @examples
#' raw_stages <- c("bud", "Bud", "floral bud", "Stage 1",
#'                 "open flower", "OF", "full bloom", "senescent")
#' standardizeStage(raw_stages)
standardizeStage <- function(stage_vector) {

  stage_lower <- tolower(trimws(stage_vector))
  result      <- stage_vector

  # Controlled vocabulary mapping
  stage_map <- list(
    "bud" = c("bud", "floral bud", "flower bud", "stage 1",
              "stage i", "stage1", "stagei", "closed bud",
              "green bud", "tight bud", "fb", "s1"),
    "opening" = c("opening", "partially open", "partially opened",
                  "opening bud", "stage 2", "stage ii",
                  "stage2", "stageii", "s2"),
    "open_flower" = c("open flower", "open", "fully open",
                      "full bloom", "anthesis", "stage 3",
                      "stage iii", "stage3", "stageiii",
                      "of", "s3", "bloom"),
    "mature_flower" = c("mature", "mature flower", "stage 4",
                        "stage iv", "stage4", "stageiv",
                        "s4", "late bloom"),
    "senescent" = c("senescent", "senescence", "wilting",
                    "wilted", "stage 5", "stage v",
                    "stage5", "stagev", "s5", "old flower"),
    "petal" = c("petal", "petals", "petal tissue"),
    "abscission_zone" = c("abscission zone", "abscission",
                          "az", "petal abscission zone",
                          "petal az"),
    "leaf" = c("leaf", "leaves", "foliage", "lamina"),
    "root" = c("root", "roots", "root tissue"),
    "stem" = c("stem", "stems", "internode"),
    "meristem" = c("meristem", "shoot apical meristem",
                   "sam", "apical meristem"),
    "stamen" = c("stamen", "stamens", "anther", "anthers",
                 "filament"),
    "pistil" = c("pistil", "pistils", "carpel", "ovary",
                 "pistil ovary")
  )

  for (std_name in names(stage_map)) {
    patterns <- stage_map[[std_name]]
    matches  <- stage_lower %in% patterns
    result[matches] <- std_name
  }

  n_standardized <- sum(result != stage_vector, na.rm = TRUE)
  message("[OrnTx] Standardized ", n_standardized,
          " of ", length(stage_vector), " stage names")

  result
}
