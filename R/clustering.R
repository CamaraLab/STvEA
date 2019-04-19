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


#' Perform louvain clustering on CODEX data
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


#' Run HDBSCAN over parameter ranges
#' Must use Python HDBSCAN because it has 2 important parameters
#' (min_cluster_size and min_samples) but R HDBSCAN only has one (minPts)
#'
#' @export
#'
ParameterScan <- function(stvea_object, min_cluster_size_range, min_sample_range,
                          ...,
                          python_dir="/usr/local/lib/R/site-library/STvEA/python",
                          cache_dir=NULL) {
  if (is.null(stvea_object@cite_latent)) {
    stop("stvea_object must contain CITE-seq latent space")
  }
  stvea_object@hdbscan_param_scan <- ParameterScan.internal(stvea_object@cite_latent,
                         min_cluster_size_range, min_sample_range,
                         ...,
                         python_dir="/usr/local/lib/R/site-library/STvEA/python",
                         cache_dir=NULL)
  return(stvea_object)
}


#' Keep only HDBSCAN results that pass a certain silhouette score cutoff
#' and create a dissimilarity matrix between cells from the clusterings
#' Perform agglomerative hierarchical clustering on the
#' consensus dissimilarity matrix
#'
#' @export
#'
ConsensusCluster <- function(stvea_object, silhouette_cutoff, num_cluster) {
  if (is.null(stvea_object@hdbscan_param_scan)) {
    stop("Please run ParameterScan first")
  }
  stvea_object@cite_clusters <- ConsensusCluster.internal(stvea_object@hdbscan_param_scan,
                                                         stvea_object@cite_latent,
                                                         silhouette_cutoff,
                                                         num_cluster)
  return(stvea_object)
}


# Functions with matrix parameters, not using STvEA.data class

#' Perform louvain clustering on CODEX data
#' Takes matrices and data frames instead of STvEA.data class
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


#' Run HDBSCAN over parameter ranges
#' Must use Python HDBSCAN because it has 2 important parameters
#' (min_cluster_size and min_samples) but R HDBSCAN only has one (minPts)
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @import rPython
#' @importFrom umap umap
#'
#' @export
#'
ParameterScan.internal <- function(cite_latent, min_cluster_size_range, min_sample_range,
                       ...,
                       python_dir="/usr/local/lib/R/site-library/STvEA/python",
                       cache_dir=NULL) {
  cite_latent_tmp <- as.matrix(cite_latent)
  colnames(cite_latent_tmp) <- NULL
  cat("Running UMAP on the CITE-seq latent space\n")
  umap_latent <- umap(cite_latent, n_components = ncol(cite_latent), ...)$layout
  umap_latent_tmp <- as.matrix(umap_latent)
  colnames(umap_latent_tmp) <- NULL
  rPython::python.load(paste(python_dir,"/consensus_clustering.py", sep=""))
  cat("Running HDBSCAN on the UMAP space\n")
  if (is.null(cache_dir)) {
    hdbscan_labels <- rPython::python.call("run_hdbscan", cite_latent_tmp, umap_latent_tmp, min_cluster_size_range, min_sample_range)
  } else {
    hdbscan_labels <- rPython::python.call("run_hdbscan", cite_latent_tmp, umap_latent_tmp, min_cluster_size_range, min_sample_range, cache_dir)
  }
  all_scores <- NULL
  hdbscan_results <- list()
  for (i in 1:length(hdbscan_labels)) {
    score <- mean(silhouette(x=hdbscan_labels[[i]], dist = dist(cite_latent))[,3])
    all_scores <- c(all_scores, score)
    hdbscan_results[[i]] <- list("cluster_labels" = hdbscan_labels[[i]], "silhouette_score" = score)
  }
  hist(all_scores, breaks=100, main="Histogram of silhouette scores", xlab = "Silhouette score", ylab = "Number of calls to HDBSCAN")
  return(hdbscan_results)
}


#' Keep only HDBSCAN results that pass a certain silhouette score cutoff
#' and create a dissimilarity matrix between cells from the clusterings
#' Perform agglomerative hierarchical clustering on the
#' consensus dissimilarity matrix
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @export
#'
ConsensusCluster.internal <- function(hdbscan_results, cite_latent, silhouette_cutoff, num_cluster) {
  num_cells <- length(hdbscan_results[[1]]$cluster_labels)
  consensus_matrix <- matrix(0, nrow=num_cells, ncol=num_cells)
  total_runs <- 0

  for (result in hdbscan_results) {
    if (result$silhouette_score >= silhouette_cutoff) {
      sim_matrix <- matrix(0, nrow=num_cells, ncol=num_cells)
      for (cell1 in 1:num_cells) {
        if (result$cluster_labels[cell1] != -1) {
          sim_matrix[,cell1] <- 1*(result$cluster_labels == result$cluster_labels[cell1])
        }
      }
      diag(sim_matrix) <- rep(1,num_cells)
      consensus_matrix <- consensus_matrix - sim_matrix
      total_runs <- total_runs + 1
    }
  }

  consensus_matrix <- consensus_matrix + total_runs
  if (total_runs > 0) {
    consensus_matrix <- consensus_matrix / total_runs
  }

  hier_tree <- hclust(as.dist(consensus_matrix), method="average")
  plot(hier_tree)
  consensus_labels <- cutree(hier_tree, k=num_cluster)
  return(consensus_labels)
}


#' Assign user input labels to clusters for future plotting
#'
clustering_names <- function(labels) {
}
