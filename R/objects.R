
#' Class to hold all data matrices needed for mapping CODEX to CITE-seq
#'
MappingDataHolder <- setClass(
  Class = 'MappingDataHolder',
  slots = c(
    cite_latent = 'ANY', # cite x low
    cite_protein = 'ANY', # cite x 30
    cite_clean = 'ANY', # cite x 30
    cite_norm = 'ANY', #cite x 30

    codex_protein = 'ANY', # codex x 30
    codex_size = 'numeric',
    codex_blanks = 'ANY',
    codex_filter = 'ANY', # codex x 30
    codex_clean = 'ANY', # codex x 30

    corrected_codex = 'ANY', # codex x 30 (can only do this if CITE always ref)
    codex_nn = 'dgCMatrix', # cite x codex
    cite_nn = 'dgCMatrix' # codex x cite
  )
)

#' Check input data and create MappingDataHolder
#'
#' @param cite_protein Raw expression data for CITE-seq proteins (cell x protein)
#' @param cite_latent Low dimensional latent space for CITE-seq mRNA (cell x dim)
#' @param codex_protein Segmented and spillover corrected CODEX protein expression (cell x protein)
#' @param codex_size Size of each CODEX cell (vector)
#' @param codex_blanks Segmented and spillover corrected CODEX blank channel expression (cell x blank channel)
#'
#' @export
#'
SetDataMapping <- function(cite_protein, cite_latent, codex_protein, codex_size, codex_blanks) {
  if (any(colnames(cite_protein) != colnames(codex_protein))) {
    stop("CITE-seq and CODEX datasets must have the same proteins in the same order")
  }
  if (nrow(cite_protein) != nrow(cite_latent)) {
    stop("CITE-seq protein and latent space must have same number of cells")
  }
  if (nrow(codex_protein) != length(codex_size) || nrow(codex_blanks) != length(codex_size)) {
    stop("CODEX protein, size, and blanks must have same number of cells")
  }
  mapping_object <- new(
    Class = "MappingDataHolder",
    cite_protein = cite_protein,
    cite_latent = cite_latent,
    codex_protein = codex_protein,
    codex_size = codex_size,
    codex_blanks = codex_blanks
  )
  return(mapping_object)
}


# other matrices: cite_gene, codex_spatial, codex_gene,
# cite_cluster_labels, codex_cluster_labels, cite_umap_embedding

# maybe have easy method for transfering data between holders?
