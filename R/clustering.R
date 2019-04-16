
# Function wrappers using STvEA.data class

#' Perform HDBSCAN consensus clustering on CITE-seq latent space
#'
#' @export
#'
ClusterCITE <- function(stvea_object) {
  if (is.null(stvea_object@cite_latent)) {
    stop("stvea_object must contain CITE-seq latent space")
  }
}


#' Get k nearest neighbors for CODEX protein space. If CleanCODEX
#' has been run, uses cleaned CODEX protein data, otherwise uses
#' original CODEX protein data.
#'
#' @param stvea_object STvEA.data class object with CODEX protein data
#' @param k number of nearest neighbors to find
#'
#' @export
#'
KnnCODEX <- function(stvea_object, k=5) {
  if (!is.null(stvea_object$codex_clean)) {
    stvea_object@codex_knn <- CorKNN(stvea_object$codex_clean, k=k)
  } else if (!is.null(stvea_object$codex_protein)) {
    stvea_object@codex_knn <- CorKNN(stvea_object$codex_protein, k=k)
  } else {
    stop("stvea_object must contain either cleaned or original CODEX protein data")
  }
  return(stvea_object)
}


#' Calls ClusterCODEX.internal with STvEA.data object
#'
#' @param stvea_object STvEA.data class object with KNN matrix for CODEX protein space
#' @param k number of nearest neighbors to use in graph-based clustering
#' If NULL, will use all nearest neighbors from the KNN matrix.
#'
#' @export
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


# Functions with matrix parameters, not using STvEA.data class

#' Perform louvain clustering on CODEX data
#'
#' @param codex_knn matrix of (codex cells x k) nearest neighbor indices
#' in the CODEX protein space
#' @param k number of nearest neighbors to use in graph-based clustering
#'
#' @importFrom igraph multilevel.community graph_from_edgelist
#' @import Matrix
#'
#' @export
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


#' Perform HDBSCAN consensus clustering on CITE-seq latent space
#'
ClusterCITE.internal <- function(cite_latent) {
  umap_out <- umap(cite_latent, n_components = ncol(cite_latent), ...)$layout

}


#' Run HDBSCAN over parameter ranges
#' Must use Python HDBSCAN because it has 2 important parameters
#' (min_cluster_size and min_samples) but R HDBSCAN only has one (minPts)
#'
#' @import rPython
#'
#' @export
#'
RunHDBSCAN <- function(cite_latent, umap_latent, min_cluster_size_range, min_sample_range, cache_dir=NULL) {
  cite_latent_tmp <- as.matrix(cite_latent)
  colnames(cite_latent_tmp) <- NULL
  #umap_latent <- umap(cite_latent, n_components = ncol(cite_latent), ...)$layout
  umap_latent_tmp <- as.matrix(umap_latent)
  colnames(umap_latent_tmp) <- NULL
  rPython::python.load('inst/python/consensus_clustering.py')
  if (is.null(cache_dir)) {
    rPython::python.call("run_hdbscan", cite_latent_tmp, cite_latent_tmp, min_cluster_size_range, min_sample_range)
  } else {
    rPython::python.call("run_hdbscan", cite_latent_tmp, cite_latent_tmp, min_cluster_size_range, min_sample_range, cache_dir)
  }
}


#' Keep only HDBSCAN results that pass a certain silhouette score cutoff
#' and create a dissimilarity matrix between cells from the clusterings
#'
#' @export
#'
GetConsensusMatrix <- function(hdbscan_results, silhouette_cutoff) {
  num_cells <- hdbscan_results[[1]]$num_cells
  dissim_matrix <- Matrix(0, nrow = num_cells, ncol = num_cells, sparse = TRUE)
  all_scores <- NULL
  total_runs <- 0
  for (result in hdbscan_results) {
    if (result$score >= silhouette_cutoff) {
      result_matrix <- sparseMatrix(i = result$rows+1, j = result$cols+1,
                                   x = result$data, dims = c(num_cells, num_cells))
      dissim_matrix <- dissim_matrix + result_matrix
      all_scores <- c(all_scores, result$score)
      total_runs <- total_runs + 1
    }
  }
  dissim_matrix <- as.matrix(dissim_matrix)
  dissim_matrix <- dissim_matrix + total_runs
  if (total_runs > 0) {
    dissim_matrix <- dissim_matrix / total_runs
  }
  return(list("dissim_matrix" = dissim_matrix, "all_scores" = all_scores))
}


#' Perform agglomerative hierarchical clustering on the
#' dissimilarity matrix output by GetConsensusMatrix
#'
#' @export
#'
ConsensusCluster <- function(dissim_matrix, num_cluster) {
  hier_tree <- hclust(dissim_matrix, method="average")
  clust_labels <- cutree(hier_tree, k=num_cluster)
  return(clust_labels)
}


#' Assign user input labels to clusters for future plotting
#'
clustering_names <- function(labels) {
}
