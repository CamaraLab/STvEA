
#' Class to hold all data matrices needed for mapping CODEX to CITE-seq
#' and following analysis
#'
STvEA.data <- setClass(
  Class = 'STvEA.data',
  slots = c(
    cite_mRNA = 'ANY', # cite_cells x genes
    cite_latent = 'ANY', # cite cells x low
    cite_emb = 'ANY', # cite cells x 2
    cite_clusters = 'vector',
    cite_protein = 'ANY', # cite cells x proteins
    cite_clean = 'ANY', # cite cells x proteins
    cite_norm = 'ANY', # cite cells x proteins

    codex_protein = 'ANY', # codex cells x proteins
    codex_size = 'numeric',
    codex_blanks = 'ANY', # codex cells x blank channels
    codex_spatial = 'ANY', # codex cells x 3
    codex_clusters = 'vector',
    codex_emb = 'ANY', # codex cells x 2
    codex_clean = 'ANY', # codex cells x proteins

    corrected_codex = 'ANY', # codex cells x proteins
    codex_nn = 'dgCMatrix', # cite cells x codex cells
    cite_nn = 'dgCMatrix', # codex cells x cite cells
    codex_mRNA = 'ANY' # codex cells x genes
  )
)


#' Check input data and create STvEA.data object
#'
#' @param cite_protein Raw expression data for CITE-seq proteins (cell x protein)
#' @param cite_latent Low dimensional latent space for CITE-seq mRNA (cell x dim)
#' @param codex_protein Segmented and spillover corrected CODEX protein expression (cell x protein)
#' @param codex_size Size of each CODEX cell (vector)
#' @param codex_blanks Segmented and spillover corrected CODEX blank channel expression (cell x blank channel)
#'
#' @export
#'
SetData <- function(cite_protein, cite_latent, codex_protein, codex_size, codex_blanks) {
  if (any(colnames(cite_protein) != colnames(codex_protein))) {
    stop("CITE-seq and CODEX datasets must have the same proteins in the same order")
  }
  if (nrow(cite_protein) != nrow(cite_latent)) {
    stop("CITE-seq protein and latent space must have same number of cells")
  }
  if (nrow(codex_protein) != length(codex_size) || nrow(codex_blanks) != length(codex_size)) {
    stop("CODEX protein, size, and blanks must have same number of cells")
  }
  stvea_object <- new(
    Class = "STvEA.data",
    cite_protein = cite_protein,
    cite_latent = cite_latent,
    codex_protein = codex_protein,
    codex_size = codex_size,
    codex_blanks = codex_blanks
  )
  return(stvea_object)
}


# other matrices: cite_gene, codex_spatial, codex_gene,
# cite_cluster_labels, codex_cluster_labels, cite_umap_embedding

# maybe have easy method for transfering data between holders?
