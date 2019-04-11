#' This file contains functions taken directly from the Seurat V3 source code,
#' edited to take matrices as input rather than Seurat objects and to be run on
#' only two datasets
#'
#' Modified from Seurat package
#' @references Stuart et al, bioRxiv 2018
#' @references \url{https://github.com/satijalab/seurat/tree/release/3.0}


#' Uses CCA to merge two datasets into a shared space
#'
#' @param object1 a (k features x n cells) matrix
#' @param object2 a (k features x m cells) matrix
#' @param standardize Standardize matrices - scales columns to have unit variance
#' and mean 0
#' @param num.cc Number of canonical vectors to calculate
#' @param verbose ...
#'
#' @return a ((n+m) cells x num.cc features) matrix
#'
#' @importFrom irlba irlba
#'
#' @export
#'
RunCCA <- function(
  object1,
  object2,
  standardize = FALSE,
  l2.norm = FALSE,
  num.cc = 20,
  setseed = TRUE,
  ...
) {
  if (setseed) {
    set.seed(seed = 42)
  }
  cells1 <- colnames(x = object1)
  cells2 <- colnames(x = object2)
  if (standardize) {
    object1 <- scale(x = object1, TRUE, apply(object1,2,sd))
    object2 <- scale(x = object2, TRUE, apply(object2,2,sd))
  }
  mat3 <- crossprod(x = object1, y = object2)

  cca.svd <- irlba(A = mat3, nv = num.cc)
  cca.data <- rbind(cca.svd$u, cca.svd$v)
  colnames(x = cca.data) <- paste0("CC", 1:num.cc)
  rownames(cca.data) <- c(cells1, cells2)
  cca.data <- apply(
    X = cca.data,
    MARGIN = 2,
    FUN = function(x) {
      if (sign(x[1]) == -1) {
        x <- x * -1
      }
      return(x)
    }
  )
  if (l2.norm){
    cca.data <- L2Norm(mat = cca.data)
  }
  return(list(ccv = cca.data, d = cca.svd$d))
}


#' Find nearest neighbors
#'
#' @param ref_emb a (cell x feature) embedding of a protein expression matrix
#' @param query_emb a (cell x feature) embedding of protein expression matrix to be corrected
#' @param k number of nearest neighbor to find between the matrices
#' @param eps error bound on nearest neighbor search (see eps parameter for nn2)
#' @param verbose ...
#'
#' @return list of nearest neighbor lists:
#' nnab = for each cell in a, its knn in b
#' nnba = for each cell in b, its knn in a
#' nnbb = knn within b, nnaa = knn within a
#'
#' @importFrom RANN nn2
#'
#' @export
#'
FindNN <- function(
  ref_emb,
  query_emb,
  k = 300,
  eps = 0,
  verbose = FALSE
) {
  if (verbose) {
    message("Finding neighborhoods")
  }
  nn.rr <- nn2(
    data = ref_emb,
    k = k + 1,
    eps = eps
  )
  nn.qq <- nn2(
    data = query_emb,
    k = k + 1,
    eps = eps
  )
  nn.rq <- nn2(
    data = query_emb,
    query = ref_emb,
    k = k,
    eps = eps
  )
  nn.qr <- nn2(
    data = ref_emb,
    query = query_emb,
    k = k,
    eps = eps
  )
  return(list('nn.rr' = nn.rr, 'nn.rq' = nn.rq, 'nn.qr' = nn.qr, 'nn.qq' = nn.qq,
              'cellsr' = row.names(ref_emb), 'cellsq' = row.names(query_emb)))
}


#' Find nearest neighbors
#'
#' @param ref_emb a (cell x feature) embedding of a protein expression matrix
#' @param query_emb a (cell x feature) embedding of protein expression matrix to be corrected
#' @param rna_mat a (cell x feature) embedding of the mRNA expression matrix from CITE-seq
#' @param cite_index which matrix (1 or 2) is the CITE-seq protein expression matrix
#' @param k number of nearest neighbor to find between the matrices
#' @param eps error bound on nearest neighbor search (see eps parameter for nn2)
#' @param verbose ...
#'
#' @return list of nearest neighbor lists:
#' nnab = for each cell in a, its knn in b
#' nnba = for each cell in b, its knn in a
#' nnbb = knn within b, nnaa = knn within a
#'
#' @importFrom RANN nn2
#'
#' @export
#'
FindNNrna <- function(
  ref_emb,
  query_emb,
  rna_mat,
  cite_index = 1,
  k = 300,
  eps = 0,
  verbose = FALSE
) {
  if (verbose) {
    message("Finding neighborhoods")
  }

  if (cite_index == 1) {
    nn.rr <- CorNN(
      data = rna_mat,
      k = k + 1
    )
    nn.qq <- nn2(
      data = query_emb,
      k = k + 1,
      eps = eps
    )
  } else {
    nn.rr <- nn2(
      data = ref_emb,
      k = k + 1,
      eps = eps
    )
    nn.qq <- CorNN(
      data = rna_mat,
      k = k + 1
    )
  }
  nn.rq <- nn2(
    data = query_emb,
    query = ref_emb,
    k = k,
    eps = eps
  )
  nn.qr <- nn2(
    data = ref_emb,
    query = query_emb,
    k = k,
    eps = eps
  )
  return(list('nn.rr' = nn.rr, 'nn.rq' = nn.rq, 'nn.qr' = nn.qr, 'nn.qq' = nn.qq,
              'cellsr' = row.names(ref_emb), 'cellsq' = row.names(query_emb)))
}


#' Find Anchor pairs as MNN between two datasets
#'
#' @param neighbors a list of knn neighbors between two datasets as output by FindNN
#' @param k.anchor number of nearest neighbors to look through for anchors
#' @param verbose ...
#'
#' @return a matrix where each row is an anchor pair of cell indices
#'
#' @export
#'
FindAnchorPairs <- function(
  neighbors,
  k.anchor = 5,
  verbose = FALSE
) {
  if (verbose) {
    message("Finding mutual nearest neighborhoods")
  }

  max.nn <- c(ncol(x = neighbors$nn.rq$nn.idx), ncol(x = neighbors$nn.qr$nn.idx))
  if (any(k.anchor > max.nn)) {
    warning('Requested k.anchor = ', k.anchor, ', only ', min(max.nn), ' in dataset')
    k.anchor <- min(max.nn)
  }

  ncell <- 1:nrow(x = neighbors$nn.rq$nn.idx)
  anchors <- list()
  # pre allocate vector
  anchors$cellr <- rep(x = 0, length(x = ncell) * 5)
  anchors$cellq <- anchors$cellr
  anchors$score <- anchors$cellr + 1
  idx <- 0
  for (cell in ncell) {
    neighbors.rq <- neighbors$nn.rq$nn.idx[cell, 1:k.anchor]
    mutual.neighbors <- which(
      x = neighbors$nn.qr$nn.idx[neighbors.rq, 1:k.anchor, drop = FALSE] == cell,
      arr.ind = TRUE
    )[, 1]
    for (i in neighbors.rq[mutual.neighbors]){
      idx <- idx + 1
      anchors$cellr[idx] <- cell
      anchors$cellq[idx] <- i
      anchors$score[idx] <- 1
    }
  }
  anchors$cellr <- anchors$cellr[1:idx]
  anchors$cellq <- anchors$cellq[1:idx]
  anchors$score <- anchors$score[1:idx]
  anchors <- t(x = do.call(what = rbind, args = anchors))
  anchors <- as.matrix(x = anchors)
  return(anchors)
}


#' Filter anchors, keeping those that are also in the knn
#'  of the L2Norm of the original feature space
#'
#' @param ref_mat a (cell x feature) protein expression matrix
#' @param query_mat a (cell x feature) protein expression matrix to be corrected
#' @param anchors the list of anchors between matrix1 and matrix2 as output by FindAnchorPairs
#' @param k.filter number of nearest neighbors in original feature space to use for filtering
#' @param eps error bound on the neighbor finding algorithm (from \code{\link{RANN}})
#' @param verbose ...
#'
#' @return a matrix where each row is an anchor pair of cell indices
#'
#' @export
#'
FilterAnchors <- function(
  ref_mat,
  query_mat,
  anchors,
  k.filter = 200,
  eps = 0,
  verbose = FALSE
) {
  if (verbose) {
    message("Filtering Anchors")
  }

  nn1 <- CorNN(
    data = query_mat,
    query = ref_mat,
    k = k.filter
  )
  nn2 <- CorNN(
    data = ref_mat,
    query = query_mat,
    k = k.filter
  )

  position1 <- sapply(X = 1:nrow(x = anchors), FUN = function(x) {
    any(anchors[x, "cellq"] == nn1$nn.idx[anchors[x, "cellr"], ])
  })
  position2 <- sapply(X = 1:nrow(x = anchors), FUN = function(x) {
    any(anchors[x, "cellr"] == nn2$nn.idx[anchors[x, "cellq"], ])
  })
  anchors <- anchors[position1 | position2, ]
  if (verbose) {
    message("\tRetained ", nrow(x = anchors), " anchors")
  }
  return(anchors)
}


#' Score anchors between 0 and 1 based on shared nearest neighbors
#'
#' @param neighbors list of neighbors in CCA space as output by FindNN
#' @param anchors list of anchors as output by FindAnchorPairs or FilterAnchors
#' @param num_cells_ref total number of cells in dataset1
#' @param num_cells_query total number of cells in dataset2
#' @param k.score number of nn to use in shared nearest neighbor scoring
#' @param verbose ...
#'
#' @import Matrix
#'
#' @return a matrix where each row is an anchor pair of cell indices with its score
#'
#' @export
#'
ScoreAnchors <- function(
  neighbors,
  anchors,
  num_cells_ref,
  num_cells_query,
  k.score = 30,
  verbose = FALSE
) {
  anchor.df <- as.data.frame(anchors)
  anchor.df$cellq <- anchor.df$cellq + num_cells_ref
  if (verbose) {
    message("Scoring anchors")
  }
  max.nn <- c(ncol(x = neighbors$nn.rr$nn.idx), ncol(x = neighbors$nn.rq$nn.idx),
              ncol(x = neighbors$nn.qr$nn.idx), ncol(x = neighbors$nn.qq$nn.idx))
  if (any(k.score > max.nn)) {
    warning('Requested k.score = ', k.score, ', only ', min(max.nn), ' in dataset')
    k.score <- min(max.nn)
  }
  total.cells <- num_cells_ref + num_cells_query
  nn.m1 <- ConstructNNMat(nn.idx = neighbors$nn.rr$nn.idx[,1:k.score],
                          offset1 = 0, offset2 = 0, dims = c(total.cells, total.cells))
  nn.m2 <- ConstructNNMat(nn.idx = neighbors$nn.rq$nn.idx[,1:k.score],
                          offset1 = 0, offset2 = num_cells_ref, dims = c(total.cells, total.cells))
  nn.m3 <- ConstructNNMat(nn.idx = neighbors$nn.qr$nn.idx[,1:k.score],
                          offset1 = num_cells_ref, offset2 = 0, dims = c(total.cells, total.cells))
  nn.m4 <- ConstructNNMat(nn.idx = neighbors$nn.qq$nn.idx[,1:k.score],
                          offset1 = num_cells_ref, offset2 = num_cells_ref, dims = c(total.cells, total.cells))
  k.matrix <- nn.m1 + nn.m2 + nn.m3 + nn.m4
  anchor.only <- sparseMatrix(i = anchor.df[, 1], j = anchor.df[, 2], x = 1, dims = c(total.cells, total.cells))

  jaccard.dist <- tcrossprod(x = k.matrix)
  anchor.matrix <- jaccard.dist * anchor.only

  anchor.matrix <- as(object = anchor.matrix, Class = "dgTMatrix")
  anchor.new <- data.frame(
    'cellr' = anchor.matrix@i + 1,
    'cellq' = anchor.matrix@j + 1,
    'score' = anchor.matrix@x
  )
  anchor.new$cellq <- anchor.new$cellq - num_cells_ref
  max.score <- quantile(anchor.new$score, 0.9)
  min.score <- quantile(anchor.new$score, 0.01)
  anchor.new$score <- anchor.new$score - min.score
  anchor.new$score <- anchor.new$score / (max.score - min.score)
  anchor.new$score[anchor.new$score > 1] <-  1
  anchor.new$score[anchor.new$score < 0] <- 0
  anchor.new <- as.matrix(x = anchor.new)
  return(anchor.new)
}


#' Construct nearest neighbor matrix from nn.idx
#'
#' @param nn.idx Nearest neighbor index matrix (nn.idx from RANN)
#' @param offset1 Offsets for the first neighbor
#' @param offset2 Offsets for the second neighbor
#'
#' @return returns a sparse matrix representing the NN matrix
#'
ConstructNNMat <- function(nn.idx, offset1, offset2, dims) {
  k <- ncol(x = nn.idx)
  j <- as.numeric(x = t(x = nn.idx)) + offset2
  i <- ((1:length(x = j)) - 1) %/% k + 1 + offset1
  nn.mat <- sparseMatrix(i = i, j = j, x = 1, dims = dims)
  return(nn.mat)
}


#' L2 normalize the columns (or rows) of a given matrix
#'
#' @param mat matrix to cosine normalize
#' @param MARGIN Perform normalization over rows (1) or columns (2)
#'
#' @return returns l2-normalized matrix
#'
L2Norm <- function(mat, MARGIN = 1){
  normalized <- sweep(
    x = mat,
    MARGIN = MARGIN,
    STATS = apply(
      X = mat,
      MARGIN = MARGIN,
      FUN = function(x){
        sqrt(x = sum(x ^ 2))
      }
    ),
    FUN = "/"
  )
  normalized[!is.finite(x = normalized)] <- 0
  return(normalized)
}


#' Calculate anchor vectors between reference and query dataset
#' This is matrix B from Methods in Stuart et al.
#'
#' @param ref_mat a (cell x feature) protein expression matrix
#' @param query_mat a (cell x feature) protein expression matrix to be corrected
#' @param neighbors a list of neighbors
#' @param anchors a list of anchors (MNN) as from FindAnchorPairs
#' @param verbose ...
#'
#' @export
#'
FindIntegrationMatrix <- function(
  ref_mat,
  query_mat,
  neighbors,
  anchors,
  verbose = FALSE
) {
  nn.cellsr <- neighbors$cellsr
  nn.cellsq <- neighbors$cellsq
  if (verbose) {
    message("Finding integration vectors")
  }
  anchors.r <- nn.cellsr[anchors[, "cellr"]]
  anchors.q <- nn.cellsq[anchors[, "cellq"]]
  data.use.r <- ref_mat[anchors.r, ]
  data.use.q <- query_mat[anchors.q, ]
  integration.matrix <- as.matrix(data.use.q - data.use.r)
  row.names(integration.matrix) <- anchors.q
  return(integration.matrix)
}


#' Calculates weights of the k nearest anchors for each cell in the query dataset
#'
#' @param neighbors list of neighbors in CCA space as output by FindNN
#' @param anchors list of anchors as output by FindAnchorPairs or FilterAnchors
#' @param query_mat a (cells x features) matrix for the query dataset (protein or mRNA)
#' @param integration.matrix the output of FindIntegrationMatrix
#' @param k number of nearest anchors to use in correction
#' @param sd.weight standard deviation of the Gaussian kernel
#' @param eps Error bound on the neighbor finding algorithm (from \code{\link{RANN}})
#' @param verbose ...
#'
#' @export
#'
FindWeights <- function(
  neighbors,
  anchors,
  query_mat,
  integration.matrix,
  k.weight = 300,
  sd.weight = 1,
  eps = 0,
  verbose = TRUE
) {
  if (verbose) {
    message("Finding integration vector weights")
  }

  nn.cellsr <- neighbors$cellsr
  nn.cellsq <- neighbors$cellsq
  anchors.cellsq <- nn.cellsq[anchors[, "cellq"]]

  # k nearest anchors
  kna_query <- CorNN(
    data = query_mat[anchors.cellsq, ],
    query = query_mat,
    k = k.weight + 1
  )

  distances <- kna_query$nn.dists[, -1]
  distances <- 1 - (distances / distances[, ncol(x = distances)])
  cell.index <- kna_query$nn.idx[, -1]

  if (verbose) {
    pb <- txtProgressBar(min = 1, max = length(x = nn.cellsq), initial = 1, style = 3, file = stderr())
  }

  dist.weights <- matrix(
    data = 0,
    nrow = nrow(x = integration.matrix),
    ncol = length(x = nn.cellsq)
  )
  for (cell in 1:length(x = nn.cellsq)) {
    wt <- distances[cell, ]
    cellnames <- anchors.cellsq[cell.index[cell, ]]
    names(x = wt) <- cellnames
    for (i in cellnames) {
      anchor.index <- which(rownames(integration.matrix) == i)
      dist.weights[anchor.index, cell] <- wt[[i]]
    }
    if (verbose) setTxtProgressBar(pb, cell)
  }
  if (verbose) cat("\n")
  dist.anchor.weight <- dist.weights * anchors[, "score"]
  weights <- 1 - exp(-1 * dist.anchor.weight / (2 * (1 / sd.weight)) ^ 2)  # Gaussian kernel
  weights <- sweep(weights, 2, Matrix::colSums(weights), "/")
  return(weights)
}


#' Correct the query dataset by subtracting the weighted anchor vectors
#'
#' @param ref_mat a (cell x feature) protein expression matrix
#' @param query_mat a (cell x feature) protein expression matrix to be corrected
#' @param integration.matrix matrix of anchor vectors (output of FindIntegrationMatrix)
#' @param weights weights of the anchors for each query cell (output of FindWeights)
#' @param verbose ...
#'
#' @return a (cell x feature) corrected data matrix
#'    the original reference data appended to the corrected query data
#'
#' @export
#'
TransformDataMatrix <- function(
  ref_mat,
  query_mat,
  integration.matrix,
  weights,
  verbose = FALSE
) {
  if(verbose) {
    message("Integrating data")
  }
  bv <-  t(weights) %*% integration.matrix
  integrated <- query_mat - bv

  new.expression <- rbind(ref_mat, integrated)
  return(new.expression)
}

