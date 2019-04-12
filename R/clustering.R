#' Perform HDBSCAN consensus clustering on CITE-seq latent space
#'
ClusterCITE <- function(cite_latent) {

}


#' Get k nearest neighbors for CODEX protein space. If CleanCODEX
#' has been run, uses cleaned CODEX protein data, otherwise uses
#' original CODEX protein data.
#'
#' @param stvea_object STvEA.data class object with CODEX protein data
#' @param k number of nearest neighbors to find
#'
KnnCODEX <- function(stvea_object, k=5) {
  if (!is.null(stvea_object$codex_clean)) {
    stvea_object@codex_knn <- CorKNN(stvea_object$codex_clean, k=k)
  } else if (!is.null(stvea_object$codex_protein)) {
    stvea_object@codex_knn <- CorKNN(stvea_object$codex_protein, k=k)
  } else {
    stop("stvea_object must contain either cleaned or original CODEX protein data")
  }
}


#' Calls ClusterCODEX.internal with STvEA.data object
#'
#' @param stvea_object STvEA.data class object with KNN matrix for CODEX protein space
#' @param k number of nearest neighbors to use in graph-based clustering
#' If NULL, will use all nearest neighbors from the KNN matrix.
#'
ClusterCODEX <- function(stvea_object, k=NULL) {
  if (is.null(stvea_object@codex_knn)) {
    stop("stvea_object must have KNN matrix for CODEX. Run KnnCODEX or GetUMAP first.")
  }
  if (is.null(k)) {
    stvea_object$codex_clusters <- ClusterCODEX.internal(stvea_object@codex_knn)
  } else {
    stvea_object$codex_clusters <- ClusterCODEX.internal(stvea_object@codex_knn, k)
  }
  return(stvea_object)
}


#' Perform louvain clustering on CODEX data
#'
#' @param codex_knn matrix of (codex cells x k) nearest neighbor indices
#' in the CODEX protein space
#' @param k number of nearest neighbors to use in graph-based clustering
#'
#' @importFrom igraph multilevel.community graph_from_edgelist
#' @import Matrix
#'
ClusterCODEX.internal <- function(codex_knn, k = ncol(codex_knn)) {
  if (k > ncol(codex_knn)) {
    stop("k must be less than or equal to number of nearest neighbors in codex_knn")
  }
  t_adj_list <- as.data.frame(t(codex_knn[,1:k]))
  colnames(t_adj_list) <- 1:nrow(codex_knn)
  edge_list <- gather(t_adj_list)
  edge_list <- apply(edge_list, 2, as.numeric)
  graph <- graph_from_edgelist(edge_list, directed=FALSE)
  communities <- multilevel.community(graph)
  reduced_clus <- factor(communities$membership)
  return(reduced_clus)
}


#' Assign user input labels to clusters for future plotting
#'
clustering_names <- function(labels) {
}
