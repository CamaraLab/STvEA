

#' Get knn in to_dataset for each cell in from_dataset
#' based on the corrected data matrix from AnchorCorrection
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
get_nn_matrix <- function(from_dataset, to_dataset, k = floor(nrow(to_dataset)*0.002), c = 0.1) {
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
#' @param ref_mat a (n cell x f feature) protein expression matrix
#' @param query_mat a (m cell x f feature) protein expression matrix to be corrected
#' @param rna_mat a (cell x feature) embedding of the mRNA expression matrix from CITE-seq
#' @param cite_index which matrix (1 or 2) is the CITE-seq protein expression matrix
#' @param num_chunks number of equal sized chunks to split CODEX dataset into for correction
#' @param num.cc number of canonical vectors to calculate
#' @param k.anchor number of nn used to find anchors via mutual nearest neighbors
#' @param k.filter number of nn in original feature space to use for filtering
#' @param k.score number of nn to use in shared nearest neighbor scoring
#' @param k.weight number of nn in original query feature space to make correction vectors
#'
#' @return a (m cell x f feature) expression matrix of the CODEX data corrected into the CITE-seq space
#'
#' @export
#'
correctChunksCODEX <- function(
  cite_protein,
  codex_protein,
  cite_latent,
  num_chunks,
  seed = NULL,
  num.cc = 29,
  k.anchor = 20,
  k.filter=100,
  k.score=80,
  k.weight=100,
  verbose=FALSE
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  random_ids <- sample(nrow(codex_protein), replace=FALSE)
  chunk_ids <- split(random_ids, cut(seq_along(random_ids), num_chunks, labels = FALSE))

  corrected_codex <- NULL
  for (i in 1:length(chunk_ids)) {
    codex_sample <- codex_protein[chunk_ids[[i]],]

    corrected_data <- AnchorCorrection(ref_mat = cite_protein,
                                       query_mat = codex_sample,
                                       rna_mat = cite_latent,
                                       cite_index = 1,
                                       num.cc = num.cc,
                                       k.anchor = k.anchor,
                                       k.filter = k.filter,
                                       k.score = k.score,
                                       k.weight = k.weight)

    corrected_codex <- rbind(corrected_codex, corrected_data[(nrow(cite_protein)+1):nrow(corrected_data),])
  }
  print(sum(row.names(corrected_codex) %in% row.names(codex_protein)))
  return(corrected_codex[row.names(codex_protein),])
}
