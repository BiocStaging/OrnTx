
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

#' Floriculture Trait Gene Database
#'
#' @description
#' A curated database of 64 genes associated with economically
#' important floriculture traits, anchored to Arabidopsis thaliana
#' gene IDs to enable cross-species queries regardless of annotation
#' quality. Covers flower color, floral scent, senescence/vase life,
#' floral development, and transcriptional regulation.
#'
#' @format A data frame with 64 rows and 9 variables:
#' \describe{
#'   \item{gene_symbol}{Gene symbol (e.g., "CHS", "DFR", "FT")}
#'   \item{arabidopsis_id}{Arabidopsis thaliana AGI locus ID}
#'   \item{trait}{Specific trait subcategory}
#'   \item{trait_category}{Broad trait category}
#'   \item{pathway}{Biochemical pathway name}
#'   \item{rose_evidence}{Evidence level in Rosa chinensis}
#'   \item{marigold_evidence}{Evidence level in Tagetes erecta}
#'   \item{key_pmid}{Key PubMed reference ID}
#' }
#'
#' @source Literature-curated from published floriculture
#'   transcriptomics studies. See key_pmid column for references.
#'
#' @name floriculture_traits
#' @keywords datasets
#' @usage data(floriculture_traits)
#' @return NULL
NULL

#' Marigold Carotenoid Gene Mapping
#'
#' @description
#' A curated mapping of Tagetes erecta gene IDs to carotenoid
#' biosynthesis gene symbols and Arabidopsis orthologs. Based on
#' differential expression analysis between orange (V-01) and
#' yellow (V-01M) flower color genotypes.
#'
#' @format A data frame with 4 rows and 3 variables:
#' \describe{
#'   \item{gene_id}{Tagetes erecta gene ID (QVD17 format)}
#'   \item{gene_symbol}{Carotenoid gene symbol (CCS, BCH, LCYe, PSY)}
#'   \item{arabidopsis_id}{Arabidopsis thaliana ortholog AGI ID}
#' }
#'
#' @source Differential expression analysis of NCBI BioProject
#'   PRJNA562616 using OrnAtlas runDE().
#'
#' @name marigold_carotenoid_genes
#' @keywords datasets
#' @usage data(marigold_carotenoid_genes)
#' @return NULL
NULL
