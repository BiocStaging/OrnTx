
#' Rosa chinensis Gene Annotation Database
#'
#' @description
#' A curated annotation table mapping generic gene IDs from the
#' RchiOBHm-V2 reference genome to official Rosa chinensis locus
#' names and gene biotype classifications.
#'
#' @format A data frame with 51,302 rows and 3 variables:
#' \describe{
#'   \item{gene_id}{Generic gene identifier (e.g., "gene0", "gene1")}
#'   \item{locus_name}{Official RchiOBHm-V2 locus name}
#'   \item{gene_biotype}{Gene biotype (protein_coding, lncRNA, etc.)}
#' }
#'
#' @source NCBI RefSeq GCF_002994745.2 (RchiOBHm-V2)
#'
#' @references
#' Hibrand Saint-Oyant L, et al. (2018) Nature Plants 4: 473-484.
#'
#' @name rosa_annotation
#' @keywords datasets
#' @usage data(rosa_annotation)
#' @return NULL
NULL

#' Rosa chinensis Example Expression Dataset
#'
#' @description
#' A small subset of the OrnAtlas Rosa chinensis expression atlas
#' for use in package examples and testing. Contains 500 highly
#' variable genes across 6 samples: 3 petal and 3 abscission zone.
#'
#' @format A SummarizedExperiment with 500 genes x 6 samples.
#' \describe{
#'   \item{assays}{counts: raw read counts matrix}
#'   \item{colData}{sample metadata including tissue_simple}
#' }
#'
#' @source NCBI SRA BioProjects PRJNA562083 and PRJNA594099
#'
#' @name rosa_example
#' @keywords datasets
#' @usage data(rosa_example)
#' @return NULL
NULL

