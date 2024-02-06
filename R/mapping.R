# Function wrappers using STvEA data class

#' Wraps MapCODEXtoCITE.internal using STvEA.data class
#'
#' @param stvea_object STvEA.data class with CODEX and CITE-seq data
#' @param num_chunks number of equal sized chunks to split CODEX dataset into for correction
#' @param seed set.seed before randomly sampling chunks of CODEX dataset
#' @param num_cores number of cores to use in parallelized correction of CODEX dataset
#' @param num.cc number of canonical vectors to calculate. Defaults to number of proteins
#' in common - 1
#' @param k.anchor number of nn used to find anchors via mutual nearest neighbors
#' @param k.filter number of nn in original feature space to use for filtering
#' @param k.score number of nn to use in shared nearest neighbor scoring
#' @param k.weight number of nn in original query feature space to make correction vectors
#'
#' @export
#'
MapCODEXtoCITE <- function(stvea_object,
                               num_chunks,
                               seed = NULL,
                               num_cores = 1,
                               num.cc = NULL,
                               k.anchor = 20,
                               k.filter=100,
                               k.score=80,
                               k.weight=100) {
  if (is.null(stvea_object@cite_clean) || is.null(stvea_object@codex_clean)) {
    stop("NormalizeCITE and CleanCODEX must have been run on input object first", call. =FALSE)
  }
  if (is.null(stvea_object@cite_latent)) {
    stop("Input object must contain lower dimensional latent space of the CITE-seq mRNA data", call. =FALSE)
  }

  common_proteins <- colnames(stvea_object@cite_clean)[colnames(stvea_object@cite_clean) %in% colnames(stvea_object@codex_clean)]
  if (length(common_proteins) < 2) {
    stop("There are too few proteins (< 2) in common between CITE-seq and CODEX protein matrices. It is advised to have at least 10 proteins in common.")
  }
  cite_clean_subset <- stvea_object@cite_clean[,common_proteins]
  codex_clean_subset <- stvea_object@codex_clean[,common_proteins]
  if (is.null(num.cc)) {
    num.cc <- length(common_proteins)-1
  }

  stvea_object@corrected_codex <- MapCODEXtoCITE.internal(cite_clean_subset,
                                                                codex_clean_subset,
                                                                stvea_object@cite_latent,
                                                                num_chunks=num_chunks, seed=seed,
                                                                num_cores=num_cores, num.cc=num.cc,
                                                                k.anchor=k.anchor, k.filter=k.filter,
                                                                k.score=k.score, k.weight=k.weight)
  return(stvea_object)
}


#' Gets NN matrix for both CITE->CODEX and CODEX->CITE, stores in STvEA.data
#'
#' @param stvea_object STvEA.data class with mapping correction run
#' @param k.cite number of nearest CODEX neighbors to find for each CITE-seq cell
#' @param c.cite width of Gaussian kernel in CITE-seq dataset
#' @param transfer_rna if set to TRUE (default), also transfers the mRNA data from CITE-seq to CODEX. Requires extra memory.
#'
#' @export
#'
GetTransferMatrix <- function(stvea_object,
                        k.cite = floor(nrow(stvea_object@corrected_codex)*0.002),
                        c.cite = 0.1,
                        transfer.rna = TRUE) {
  if (is.null(stvea_object@corrected_codex)) {
    stop("MapCODEXtoCITE must be run on the input object first", call. =FALSE)
  }
  stvea_object@transfer_matrix <- GetTransferMatrix.internal(stvea_object@cite_clean[,colnames(stvea_object@corrected_codex)],
                                                 stvea_object@corrected_codex,
                                                 k=k.cite,
                                                 c=c.cite)
  if (transfer_rna) {
      stvea_object@codex_mRNA <- as.matrix(stvea_object@transfer_matrix) %*%
        as.matrix(stvea_object@cite_mRNA/rowSums(stvea_object@cite_mRNA))
  }
  if (!is.null(stvea_object@cite_clusters)) {
    cluster_ids <- as.character(unique(stvea_object@cite_clusters))
    cluster_ids <- cluster_ids[order(cluster_ids)]
    cluster_matrix <- t(sapply(cluster_ids, function(x) (stvea_object@cite_clusters==x)*1))
    row.names(cluster_matrix) <- cluster_ids
    codex_cluster_matrix <- stvea_object@transfer_matrix %*% t(cluster_matrix)
    stvea_object@codex_clusters <- c(-1, cluster_ids)[(apply(codex_cluster_matrix, 1, max)!=0)*apply(codex_cluster_matrix, 1, which.max)+1]
  }
  return(stvea_object)
}



# Functions with matrix parameters, not using STvEA.data object

#' Find the knn indices in data for all cells in query
#' using Pearson's correlation dissimilarity
#'
#' Note: "query" here is used to mimic RANN documentation, not Seurat documentation
#'
#' @param data a (M cell x d feature) matrix
#' @param query a (N cell x d feature) matrix. If not specified, query = data.
#' @param k number of nearest neighbors to find
#'
#' @return a N x k matrix containing nearest neighbors for each cell in query
#'
CorNN <- function(
  data,
  query = data,
  k = 5
) {
  t_data <- t(data)
  query <- as.matrix(query)
  neighbors <- matrix(rep(0, k*nrow(query)), ncol=k)
  distances <- matrix(rep(0, k*nrow(query)), ncol=k)
  for (i in 1:nrow(query)) {
    cor_dist <- 1-cor(query[i,], t_data)
    idx <- order(cor_dist)[1:k]
    neighbors[i,] <- idx
    distances[i,] <- cor_dist[idx]
  }
  return(list(nn.idx=neighbors, nn.dists=distances))
}


#' Run all methods for anchor correction
#'
#' @param ref_mat a (cell x feature) protein expression matrix
#' @param query_mat a (cell x feature) protein expression matrix to be corrected
#' @param rna_mat a (cell x feature) embedding of the mRNA expression matrix from CITE-seq
#' @param cite_index which matrix (1 or 2) is the CITE-seq protein expression matrix
#' @param num.cc number of canonical vectors to calculate. Defaults to number of proteins - 1
#' @param k.anchor number of nn used to find anchors via mutual nearest neighbors
#' @param k.filter number of nn in original feature space to use for filtering
#' @param k.score number of nn to use in shared nearest neighbor scoring
#' @param k.weight number of nn in original query feature space to make correction vectors
#' @param verbose print status at each step
#'
#' @export
#'
AnchorCorrection <- function(
  ref_mat,
  query_mat,
  rna_mat,
  cite_index,
  num.cc = NULL,
  k.anchor = 20,
  k.filter=100,
  k.score=80,
  k.weight=100,
  verbose=FALSE
) {
  if (is.null(num.cc)) {
    num.cc <- ncol(ref_mat)-1
  }
  cca_matrix <- RunCCA(t(ref_mat), t(query_mat), standardize=TRUE, num.cc=num.cc)$ccv
  neighbors <- FindNNrna(ref_emb = cca_matrix[1:nrow(ref_mat),],
                         query_emb = cca_matrix[(nrow(ref_mat)+1):nrow(cca_matrix),],
                         rna_mat = rna_mat,
                         cite_index = cite_index,
                         k=max(k.anchor, k.score), verbose=verbose)
  anchors <- FindAnchorPairs(neighbors, k.anchor=k.anchor)
  filteredAnchors <- FilterAnchors(ref_mat, query_mat, anchors, k.filter=k.filter, verbose=verbose)
  scoredAnchors <- ScoreAnchors(neighbors, filteredAnchors, nrow(ref_mat), nrow(query_mat), k.score=k.score, verbose=verbose)
  integration.matrix <- FindIntegrationMatrix(ref_mat, query_mat, neighbors, scoredAnchors, verbose=verbose)
  weights <- FindWeights(neighbors, scoredAnchors, query_mat, integration.matrix, k.weight=k.weight, verbose=verbose)
  corrected_data <- TransformDataMatrix(ref_mat, query_mat, integration.matrix, weights, verbose=verbose)
  return(corrected_data)
}


#' Run all methods for anchor correction to map the CODEX dataset to the CITE-seq dataset
#' Since the CODEX dataset is usually much bigger, allows to subsample equal-sized
#' chunks of the CODEX dataset to map individually
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param cite_protein a (n cell x f feature) protein expression matrix
#' @param codex_protein a (m cell x f feature) protein expression matrix to be corrected
#' @param cite_latent a (cell x feature) embedding of the mRNA expression matrix from CITE-seq
#' @param num_chunks number of equal sized chunks to split CODEX dataset into for correction
#' @param seed set.seed before randomly sampling chunks of CODEX dataset
#' @param num_cores number of cores to use in parallelized correction of CODEX dataset.
#' On Windows, this must be set to 1.
#' @param num.cc number of canonical vectors to calculate. Defaults to number of proteins - 1
#' @param k.anchor number of nn used to find anchors via mutual nearest neighbors
#' @param k.filter number of nn in original feature space to use for filtering
#' @param k.score number of nn to use in shared nearest neighbor scoring
#' @param k.weight number of nn in original query feature space to make correction vectors
#'
#' @return a (m cell x f feature) expression matrix of the CODEX data corrected into the CITE-seq space
#'
#' @export
#'
MapCODEXtoCITE.internal <- function(
  cite_protein,
  codex_protein,
  cite_latent,
  num_chunks,
  seed = NULL,
  num_cores = 1,
  num.cc = NULL,
  k.anchor = 20,
  k.filter=100,
  k.score=80,
  k.weight=100
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(num.cc)) {
    num.cc <- ncol(cite_protein)-1
  }

  if (num_chunks == 1) {
    corrected_data <- AnchorCorrection(ref_mat = cite_protein,
                                       query_mat = codex_protein,
                                       rna_mat = cite_latent, cite_index = 1,
                                       num.cc = num.cc, k.anchor = k.anchor,
                                       k.filter = k.filter, k.score = k.score,
                                       k.weight = k.weight)
    return(corrected_data[(nrow(cite_protein)+1):nrow(corrected_data),])
  }

  random_ids <- sample(nrow(codex_protein), replace=FALSE)
  chunk_ids <- split(random_ids, cut(seq_along(random_ids), num_chunks, labels = FALSE))
  if (num_cores > 1) {
    chunk_ids <- split(random_ids, cut(seq_along(random_ids), num_chunks, labels = FALSE))
    corrected_data <- mclapply(1:length(chunk_ids),
                             function(i) AnchorCorrection(ref_mat = cite_protein,
                                                          query_mat = codex_protein[chunk_ids[[i]],],
                                                          rna_mat = cite_latent, cite_index = 1,
                                                          num.cc = num.cc, k.anchor = k.anchor,
                                                          k.filter = k.filter, k.score = k.score,
                                                          k.weight = k.weight),
                             mc.cores=num_cores)
  } else {
    corrected_data <- lapply(1:length(chunk_ids),
                               function(i) AnchorCorrection(ref_mat = cite_protein,
                                                            query_mat = codex_protein[chunk_ids[[i]],],
                                                            rna_mat = cite_latent, cite_index = 1,
                                                            num.cc = num.cc, k.anchor = k.anchor,
                                                            k.filter = k.filter, k.score = k.score,
                                                            k.weight = k.weight))
  }
  corrected_codex <- NULL
  for (i in 1:length(chunk_ids)) {
    corrected_codex <- rbind(corrected_codex,
                             corrected_data[[i]][(nrow(cite_protein)+1):nrow(corrected_data[[i]]),])
  }
  return(corrected_codex[row.names(codex_protein),])
}


#' Get knn in to_dataset for each cell in from_dataset
#' based on the corrected data matrix from AnchorCorrection
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param from_dataset a (cell x feature) matrix or dataframe
#' @param to_dataset a (cell x feature) matrix or dataframe
#' @param k number of nearest neighbors to find
#' @param c constant controls the width of the Gaussian kernel
#'
#' @import Matrix
#' @importFrom tidyr gather
#'
#' @export
#'
GetTransferMatrix.internal <- function(from_dataset, to_dataset, k = floor(nrow(to_dataset)*0.002), c = 0.1) {
  # get corrected data for each dataset (either as input or from name in object)
  # compute query knn from CorNN
  # weight each nn based on gaussian kernel of distance
  # create weighted nn matrix as sparse matrix
  # return nn matrix (maybe as part of object)
  nn_list <- CorNN(to_dataset , from_dataset, k=k)
  nn_idx <- nn_list$nn.idx
  row.names(nn_idx) <- 1:nrow(nn_idx)

  nn_dists_exp <- exp(nn_list$nn.dists/-c)
  nn_weights <- nn_dists_exp / rowSums(nn_dists_exp)
  row.names(nn_weights) <- 1:nrow(nn_idx)

  sparse_coords <- gather(as.data.frame(t(nn_idx)))
  sparse_coords$key <- as.numeric(sparse_coords$key)
  sparse_entries <- gather(as.data.frame(t(nn_weights)))
  nn_matrix <- sparseMatrix(i = sparse_coords$value,
                            j = sparse_coords$key,
                            x = sparse_entries$value,
                            dims=c(nrow(to_dataset),nrow(from_dataset)))

  colnames(nn_matrix) <- row.names(from_dataset)
  row.names(nn_matrix) <- row.names(to_dataset)
  return(nn_matrix)
}
