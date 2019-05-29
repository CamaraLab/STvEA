
#' Class to hold all data matrices needed for mapping CODEX to CITE-seq
#' and following analysis
#'
STvEA.data <- setClass(
  Class = 'STvEA.data',
  slots = c(
    cite_mRNA = 'ANY', # cite_cells x genes
    cite_mRNA_norm = 'ANY', # cite cells x genes - log(1+TPM) or seurat
    cite_latent = 'ANY', # cite cells x low
    cite_emb = 'ANY', # cite cells x 2
    hdbscan_param_scan = 'list', # output of ParameterScan
    cite_clusters = 'vector',
    cite_protein = 'ANY', # cite cells x proteins
    cite_clean = 'ANY', # cite cells x proteins
    cite_norm = 'ANY', # cite cells x proteins

    codex_protein = 'ANY', # codex cells x proteins
    codex_size = 'numeric',
    codex_blanks = 'ANY', # codex cells x blank channels
    codex_spatial = 'ANY', # codex cells x 3
    codex_emb = 'ANY', # codex cells x 2
    codex_knn = 'ANY', # codex cells x k
    codex_clusters = 'vector',
    codex_clean = 'ANY', # codex cells x proteins

    corrected_codex = 'ANY', # codex cells x proteins
    transfer_matrix = 'dgCMatrix', # codex cells x cite cells
                                 # k CODEX neighbors for each CITE-seq cell
    codex_mRNA = 'ANY' # codex cells x genes
  )
)


#' Set CODEX data in STvEA.data object
#'
#' @param codex_protein Segmented and spillover corrected CODEX protein expression (cell x protein)
#' @param codex_blanks Segmented and spillover corrected CODEX blank channel expression (cell x blank channel)
#' @param codex_size Size of each CODEX cell (vector)
#' @param codex_spatial xyz coordinates of each CODEX cell (cell x 3 coordinates)
#' @param stvea_object (optional) Pre-existing STvEA.data object to load data into.
#' If not provided, a new object is created.
#'
#' @export
#'
SetDataCODEX <- function(codex_protein,
                         codex_blanks,
                         codex_size,
                         codex_spatial,
                         stvea_object = NULL) {
  if (is.null(stvea_object)) {
    stvea_object <- new(
      Class = "STvEA.data",
      codex_protein = codex_protein,
      codex_blanks = codex_blanks,
      codex_size = codex_size,
      codex_spatial = codex_spatial
    )
  } else {
    stvea_object@codex_protein <- codex_protein
    stvea_object@codex_blanks <- codex_blanks
    stvea_object@codex_size <- codex_size
    stvea_object@codex_spatial <- codex_spatial
  }
  return(stvea_object)
}


#' Set CITE-seq data in STvEA.data object
#'
#' @param cite_mRNA Raw count data for CITE-seq mRNA (cell x gene)
#' @param cite_protein Raw expression data for CITE-seq proteins (cell x protein)
#' @param cite_mRNA_norm (optional) Normalized and scaled CITE-seq mRNA count data. (cell x gene)
#' If not provided, raw count data is divided by total per cell, then log transformed.
#' @param cite_latent (optional) Low dimensional latent space for CITE-seq mRNA (cell x dim)
#' @param stvea_object (optional) Pre-existing STvEA.data object to load data into.
#' If not provided, a new object is created.
#'
#' @export
#'
SetDataCITE <- function(cite_mRNA,
                        cite_protein,
                        cite_mRNA_norm = NULL,
                        cite_latent = NULL,
                        stvea_object = NULL) {
  if (is.null(cite_mRNA_norm)) {
    cite_mRNA_norm <- log(1 + 1000*(cite_mRNA/rowSums(cite_mRNA)))
  }
  if (is.null(stvea_object)) {
    stvea_object <- new(
      Class = "STvEA.data",
      cite_mRNA = cite_mRNA,
      cite_mRNA_norm = cite_mRNA_norm,
      cite_protein = cite_protein
    )
  } else {
    stvea_object@cite_mRNA <- cite_mRNA
    stvea_object@cite_mRNA_norm <- cite_mRNA_norm
    stvea_object@cite_protein <- cite_protein
  }
  if(!is.null(cite_latent)) {
    stvea_object@cite_latent <- cite_latent
  }
  return(stvea_object)
}


#' Transfer data from a Seurat object to STvEA.data object
#'
#' @param seurat_object Seurat class object
#' @param embedding_reduction (optional) name of the slot holding the dimension reduction
#' embedding used for visualization. Examples: "umap", "tsne"
#' @param latent_reduction (optional) name of the lot holding the dimension reduction to
#' be used for the CITE-seq mRNA latent space. Examples: "pca"
#' @param latent_dims number of dimensions to use for the latent space,
#' not used if latent_reduction is not defined
#'
#' @return STvEA.data class object
#'
#' @export
#'
TransferDataSeurat2 <- function(seurat_object,
                                embedding_reduction = NULL,
                                latent_reduction = NULL,
                                latent_dims = 20) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  stvea_object <- new(
    Class = "STvEA.data",
    cite_mRNA = t(seurat_object@raw.data),
    cite_protein = t(Seurat::GetAssayData(seurat_object, assay.type="CITE", slot="raw.data")),
    cite_clusters = Seurat::GetClusters(seurat_object)$cluster
  )
  if (!is.null(seurat_object@scale.data)) {
    stvea_object@cite_mRNA_norm <- t(seurat_object@scale.data)
  }
  if (!is.null(latent_reduction)) {
    cite_latent <- Seurat::GetDimReduction(
      object = seurat_object,
      reduction.type = latent_reduction,
      slot = 'cell.embeddings'
    )
    stvea_object@cite_latent <- cite_latent[,1:min(latent_dims,ncol(cite_latent))]
  }
  if (!is.null(embedding_reduction)) {
    stvea_object@cite_emb <- as.data.frame(Seurat::GetDimReduction(
      object = seurat_object,
      reduction.type = embedding_reduction,
      slot = 'cell.embeddings'
    ))
  }
  return(stvea_object)
}


