
#' Runs UMAP on CITE-seq latent space and CODEX protein expression
#' if they exist in the provided STvEA.data object. If CleanCODEX
#' has been run, calls UMAP on the cleaned CODEX protein expression.
#' Otherwise, runs UMAP on the original CODEX protein expression.
#'
#' @param stvea_object STvEA.data class object
#' @param ... parameters to pass into umap
#'
GetVisUMAP <- function(stvea_object, ...) {
  if (!is.null(stvea_object@cite_latent)) {
    res <- umap::umap(stvea_object@cite_latent, ...)
    stvea_object@cite_emb <- as.data.frame(res$layout)
  }
  if (!is.null(stvea_object@codex_clean)) {
    res <- umap::umap(stvea_object@codex_clean, ...)
    stvea_object@codex_emb <- as.data.frame(res$layout)
    stvea_object@codex_knn <- res$knn$indexes
  } else if (!is.null(stvea_object@codex_protein)) {
    res <- umap::umap(stvea_object@codex_protein, ...)
    stvea_object@codex_emb <- as.data.frame(res$layout)
    stvea_object@codex_knn <- res$knn$indexes
  }
  return(stvea_object)
}
