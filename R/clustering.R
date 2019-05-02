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
  if (!is.null(stvea_object@codex_clean)) {
    stvea_object@codex_knn <- CorNN(stvea_object@codex_clean, k=k)$nn.idx
  } else if (!is.null(stvea_object@codex_protein)) {
    stvea_object@codex_knn <- CorNN(stvea_object@codex_protein, k=k)$nn.idx
  } else {
    stop("stvea_object must contain either cleaned or original CODEX protein data", call. =FALSE)
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
    stop("stvea_object must have KNN matrix for CODEX. Run KnnCODEX or GetUMAP first.", call. =FALSE)
  }
  if (is.null(k)) {
    stvea_object@codex_clusters <- ClusterCODEX.internal(stvea_object@codex_knn)
  } else {
    stvea_object@codex_clusters <- ClusterCODEX.internal(stvea_object@codex_knn, k)
  }
  return(stvea_object)
}


#' Run HDBSCAN over parameter ranges
#' Must use Python HDBSCAN because it has 2 important parameters
#' (min_cluster_size and min_samples) but R HDBSCAN only has one (minPts)
#'
#' @param stvea_object STvEA.data class object containing CITE-seq latent space
#' @param min_cluster_size_range vector of min_cluster_size arguments to scan over
#' @param min_sample_range vector of min_sample arguments to scan over
#' @param ... extra parameters to be passed into UMAP
#'
#' @export
#'
ParameterScan <- function(stvea_object, min_cluster_size_range, min_sample_range, ...) {
  if (is.null(stvea_object@cite_latent)) {
    stop("stvea_object must contain CITE-seq latent space", call. =FALSE)
  }
  stvea_object@hdbscan_param_scan <- ParameterScan.internal(stvea_object@cite_latent,
                         min_cluster_size_range, min_sample_range, ...)
  return(stvea_object)
}


#' Keep only HDBSCAN results that pass a certain silhouette score cutoff
#' and create a dissimilarity matrix between cells from the clusterings
#' Perform agglomerative hierarchical clustering on the
#' consensus dissimilarity matrix
#'
#' @param stvea_object STvEA.data class object after ParameterScan has been run
#' @param silhouette_cutoff minimum silhouette score to keep clustering
#' @param inconsistent_value input parameter to python fcluster determining
#' where clusters are cut in the hierarchical tree
#' @param min_cluster_size cells in clusters smaller than this value are
#' assigned a cluster ID of -1, indicating no cluster assignment
#'
#' @export
#'
ConsensusCluster <- function(stvea_object, silhouette_cutoff, inconsistent_value, min_cluster_size) {
  if (is.null(stvea_object@hdbscan_param_scan)) {
    stop("Please run ParameterScan first", call. =FALSE)
  }
  stvea_object@cite_clusters <- ConsensusCluster.internal(stvea_object@hdbscan_param_scan,
                                                         stvea_object@cite_latent,
                                                         silhouette_cutoff,
                                                         inconsistent_value,
                                                         min_cluster_size)
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
    warning("k must be less than or equal to number of nearest neighbors from GetUmapCODEX or KnnCODEX", call. =FALSE)
    k = ncol(codex_knn)
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
#' @param cite_latent CITE-seq latent space (cells x dimensions)
#' @param min_cluster_size_range vector of min_cluster_size arguments to scan over
#' @param min_sample_range vector of min_sample arguments to scan over
#' @param ... extra parameters to be passed into UMAP
#'
#' @importFrom cluster silhouette
#'
#' @export
#'
ParameterScan.internal <- function(cite_latent, min_cluster_size_range, min_sample_range, ...) {
  if (!requireNamespace("umap", quietly = TRUE)) {
    stop("Package \"umap\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package \"reticulate\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if(!reticulate::py_module_available("numpy")) {
    stop("Python must have package \"numpy\" installed for this function to work. Please install it.",
         call.=FALSE)
  }
  if(!reticulate::py_module_available("scipy")) {
    stop("Python must have package \"scipy\" installed for this function to work. Please install it.",
         call.=FALSE)
  }
  if(!reticulate::py_module_available("hdbscan")) {
    stop("Python must have package \"hdbscan\" installed for this function to work. Please install it.",
         call.=FALSE)
  }
  if(!reticulate::py_module_available("sklearn")) {
    stop("Python must have package \"sklearn\" installed for this function to work. Please install it.",
         call.=FALSE)
  }

  cite_latent_tmp <- as.matrix(cite_latent)
  colnames(cite_latent_tmp) <- NULL
  cat("Running UMAP on the CITE-seq latent space\n")
  umap_latent <- umap::umap(cite_latent, n_components = ncol(cite_latent), ...)$layout
  umap_latent_tmp <- as.matrix(umap_latent)
  colnames(umap_latent_tmp) <- NULL

  python_dir <- paste(installed.packages()["STvEA","LibPath"],"/STvEA/python/", sep="")
  cache_dir <- paste(python_dir, "/tmp/",sep="")

  reticulate::source_python(paste(python_dir,"consensus_clustering.py",sep=""))
  cat("Running HDBSCAN on the UMAP space\n")
  hdbscan_labels <- run_hdbscan(cite_latent_tmp, umap_latent_tmp, min_cluster_size_range, min_sample_range, cache_dir)

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
#' @param hdbscan_results output of ParameterScan.internal
#' @param cite_latent CITE-seq latent space (cells x dimensions)
#' @param silhouette_cutoff minimum silhouette score to keep clustering
#' @param inconsistent_value input parameter to python fcluster determining
#' where clusters are cut in the hierarchical tree
#' @param min_cluster_size cells in clusters smaller than this value are
#' assigned a cluster ID of -1, indicating no cluster assignment
#'
#' @export
#'
ConsensusCluster.internal <- function(hdbscan_results, cite_latent, silhouette_cutoff, inconsistent_value, min_cluster_size) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package \"reticulate\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

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
  colnames(consensus_matrix) <- NULL

  reticulate::source_python(paste(installed.packages()["STvEA","LibPath"],"/STvEA/python/consensus_clustering.py",sep=""))
  consensus_clusters <- consensus_cluster(consensus_matrix, inconsistent_value, min_cluster_size)

  # relabel clusters so they are sequential, except -1 (no cluster)
  clustered_pts <- consensus_clusters != -1
  new_labels <- as.numeric(factor(consensus_clusters[clustered_pts]))
  consensus_clusters[clustered_pts] <- new_labels
  return(consensus_clusters)
}


#' Assign user input labels to clusters for future plotting
#'
clustering_names <- function(labels) {
}
