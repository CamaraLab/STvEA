#' @import AdjacencyScore
#'
NULL


#' Use the Adjacency Score to evaluate colocaliztion of all pairs
#' of CITE-seq clusters mapped to the CODEX spatial positions
#'
#' @param stvea_object STvEA.data class object containing CITE-seq
#' cluster assignments, CODEX spatial coordinates, and transfer matrix
#' from GetTransferMatrix
#' @param k number of nearest neighbors to create graph
#' from CODEX spatial information
#' @param c constant used to determine width of diffusion, must be 0 <= c
#' @param num_cores integer specifying the number of cores to be used
#' in the computation. By default only one core is used.
#' @param num_perms number of permutations used to build the null
#' distribution for each feature. By default is set to 1000.
#' @param perm_estimate boolean indicating whether Gaussian distribution
#' parameters should be determined from num_perms permutations to estimate
#' the p-value. By default is set to TRUE.
#'
#' @export
#'
AdjScoreClustersCITE <- function(
  stvea_object,
  k,
  c=0,
  num_cores=1,
  num_perms=1000,
  perm_estimate=T
) {
  knn_adj <- knn_graph(stvea_object@codex_spatial, k=k)
  AdjScoreClustersCITE.internal(stvea_object@cite_clusters,
                                stvea_object@codex_transfer,
                                knn_adj,
                                c=c,
                                num_cores=num_cores,
                                num_perms=num_perms,
                                perm_estimate=perm_estimate)
}


#' Use the Adjacency Score to evaluate colocaliztion of all pairs
#' of CODEX clusters in CODEX spatial dimensions.
#' Calls Adjacency Score with c=0 to use hypergeometric null distribution
#' speed up for mutually exclusive binary features.
#'
#' @param stvea_object STvEA.data class object containing CODEX
#' cluster assignments and CODEX spatial coordinates
#' @param k number of nearest neighbors to create graph
#' from CODEX spatial information
#' @param num_cores integer specifying the number of cores to be used
#' in the computation. By default only one core is used.
#'
#' @export
#'
AdjScoreClustersCODEX <- function(stvea_object, k, num_cores=1) {
  knn_adj <- knn_graph(stvea_object@codex_spatial, k=k)
  AdjScoreClustersCODEX.internal(stvea_object$codex_clusters,
                                 knn_adj,
                                 num_cores=num_cores)
}


#' Use the Adjacency Score to evaluate colocalization of given
#' pairs of proteins in the CODEX spatial dimensions
#'
#' @param stvea_object STvEA.data class object containing CODEX
#' protein expression and CODEX spatial coordinates
#' @param protein_pairs a 2 column matrix of protein pairs where each row
#' specifies the names of the proteins in a pair
#' @param k number of nearest neighbors to create graph
#' from CODEX spatial information
#' @param c constant used to determine width of diffusion, must be 0 <= c
#' @param num_cores integer specifying the number of cores to be used
#' in the computation. By default only one core is used.
#' @param num_perms number of permutations used to build the null
#' distribution for each feature. By default is set to 1000.
#' @param perm_estimate boolean indicating whether Gaussian distribution
#' parameters should be determined from num_perms permutations to estimate
#' the p-value. By default is set to TRUE.
#'
#' @export
#'
AdjScoreProteins <- function(
  stvea_object,
  protein_pairs,
  k,
  c=0,
  num_cores=1,
  num_perms=1000,
  perm_estimate=TRUE
) {
  knn_adj <- knn_graph(stvea_object@codex_spatial, k=k)
  AdjScoreProteins.internal(stvea_object@codex_protein,
                            protein_pairs,
                            knn_adj,
                            c=c,
                            num_cores=num_cores,
                            num_perms=num_perms,
                            perm_estimate=perm_estimate)
}


#' Use the Adjacency Score to evaluate colocalization of given
#' pairs of genes mapped to the CODEX spatial positions
#'
#' @param stvea_object STvEA.data class object with CITE-seq gene
#' expression data, CODEX spatial coordinates, and transfer matrix
#' from GetTransferMatrix
#' @param gene_pairs a 2 column matrix of gene pairs where each row
#' specifies the names of the genes in a pair
#' @param k number of nearest neighbors to create graph
#' from CODEX spatial information
#' @param c constant used to determine width of diffusion, must be 0 <= c
#' @param num_cores integer specifying the number of cores to be used
#' in the computation. By default only one core is used.
#' @param num_perms number of permutations used to build the null
#' distribution for each feature. By default is set to 1000.
#' @param perm_estimate boolean indicating whether Gaussian distribution
#' parameters should be determined from num_perms permutations to estimate
#' the p-value. By default is set to TRUE.
#'
#' @export
#'
AdjScoreGenes <- function(
  stvea_object,
  gene_pairs,
  k,
  c=0,
  num_cores=1,
  num_perms=1000,
  perm_estimate=T
) {
  knn_adj <- knn_graph(stvea_object@codex_spatial, k=k)
  AdjScoreGenes.interna(stvea_object@cite_mRNA,
                        gene_pairs,
                        knn_adj,
                        c=c,
                        num_cores=num_cores,
                        num_perms=num_perms,
                        perm_estimate=perm_estimate)
}


# Functions with matrix parameters, not using STvEA.data object

#' Use the Adjacency Score to evaluate colocaliztion of all pairs
#' of CITE-seq clusters mapped to the CODEX spatial positions
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param cite_clusters a vector of cluster IDs for the CITE-seq cells
#' @param codex_transfer a (cite-seq cells x codex cells) matrix
#' of weighted nearest neighbor assignments mapping each CODEX
#' cell to k CITE-seq cells
#' @param adj_matrix a (preferrably sparse) binary matrix of
#' adjacency between the cells in the CODEX spatial coordinates
#' @param c constant used to determine width of diffusion, must be 0 <= c
#' @param num_cores integer specifying the number of cores to be used
#' in the computation. By default only one core is used.
#' @param num_perms number of permutations used to build the null
#' distribution for each feature. By default is set to 1000.
#' @param perm_estimate boolean indicating whether Gaussian distribution
#' parameters should be determined from num_perms permutations to estimate
#' the p-value. By default is set to TRUE.
#'
#' @export
#'
AdjScoreClustersCITE.internal <- function(
  cite_clusters,
  codex_transfer,
  adj_matrix,
  c=0,
  num_cores=1,
  num_perms=1000,
  perm_estimate=T
) {
  cluster_ids <- unique(cite_clusters)
  cluster_ids <- cluster_ids[order(cluster_ids)]
  cluster_matrix <- t(sapply(cluster_ids, function(x) (cite_clusters==x)*1))
  row.names(cluster_matrix) <- cluster_ids

  cluster_pairs <- t(combn(cluster_ids,2))
  for (id in cluster_ids) {
    cluster_pairs <- rbind(cluster_pairs, c(id,id))
  }

  codex_cluster_matrix <- cluster_matrix %*% codex_transfer

  adjacency_score(adj_matrix, codex_cluster_matrix, cluster_pairs,
                  c=c, num_cores=num_cores, num_perms=num_perms,
                  perm_estimate = perm_estimate)
}


#' Use the Adjacency Score to evaluate colocaliztion of all pairs
#' of CODEX clusters in CODEX spatial dimensions.
#' Calls Adjacency Score with c=0 to use hypergeometric null distribution
#' speed up for mutually exclusive binary features.
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param codex_clusters a vector of the cluster ID for each CODEX cell
#' @param adj_matrix a (preferrably sparse) binary matrix of
#' adjacency between the cells in the CODEX spatial coordinates
#' @param num_cores integer specifying the number of cores to be used
#' in the computation. By default only one core is used.
#'
#' @export
#'
AdjScoreClustersCODEX.internal <- function(
  codex_clusters,
  adj_matrix,
  num_cores=1
) {
  cluster_ids <- unique(codex_clusters)
  cluster_ids <- cluster_ids[order(cluster_ids)]
  cluster_matrix <- t(sapply(cluster_ids, function(x) (codex_clusters==x)*1))
  row.names(cluster_matrix) <- cluster_ids

  cluster_pairs <- t(combn(cluster_ids,2))
  for (id in cluster_ids) {
    cluster_pairs <- rbind(cluster_pairs, c(id,id))
  }

  # Call adjacency score with c=0 and groupings=TRUE to use hypergeometric null distribution speed up
  adjacency_score(adj_matrix, cluster_matrix, cluster_pairs, c=0, num_cores=num_cores, groupings=TRUE)
}


#' Use the Adjacency Score to evaluate colocalization of given
#' pairs of proteins in the CODEX spatial dimensions
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param codex_protein a (cells x proteins) matrix of CODEX protein expression
#' @param protein_pairs a 2 column matrix of protein pairs where each row
#' specifies the names of the proteins in a pair
#' @param adj_matrix a (preferrably sparse) binary matrix of
#' adjacency between the cells in the CODEX spatial coordinates
#' @param c constant used to determine width of diffusion, must be 0 <= c
#' @param num_cores integer specifying the number of cores to be used
#' in the computation. By default only one core is used.
#' @param num_perms number of permutations used to build the null
#' distribution for each feature. By default is set to 1000.
#' @param perm_estimate boolean indicating whether Gaussian distribution
#' parameters should be determined from num_perms permutations to estimate
#' the p-value. By default is set to TRUE.
#'
#' @export
#'
AdjScoreProteins.internal <- function(
  codex_protein,
  protein_pairs,
  adj_matrix,
  c,
  num_cores,
  num_perms,
  perm_estimate
) {
  codex_protein_cut <- t(codex_protein[,colnames(codex_protein) %in% as.vector(protein_pairs)])
  adjacency_score(adj_matrix, codex_protein_cut, protein_pairs, c=c, num_cores=num_cores,
                  num_perms=num_perms, perm_estimate=perm_estimate)
}


#' Use the Adjacency Score to evaluate colocalization of given
#' pairs of genes mapped to the CODEX spatial positions
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param codex_mRNA a (CODEX cells x genes) matrix of the CITE-seq
#' mRNA expression mapped to the CODEX cells
#' @param gene_pairs a 2 column matrix of gene pairs where each row
#' specifies the names of the genes in a pair
#' @param adj_matrix a (preferrably sparse) binary matrix of
#' adjacency between the cells in the CODEX spatial coordinates
#' @param c constant used to determine width of diffusion, must be 0 <= c
#' @param num_cores integer specifying the number of cores to be used
#' in the computation. By default only one core is used.
#' @param num_perms number of permutations used to build the null
#' distribution for each feature. By default is set to 1000.
#' @param perm_estimate boolean indicating whether Gaussian distribution
#' parameters should be determined from num_perms permutations to estimate
#' the p-value. By default is set to TRUE.
#'
#' @export
#'
AdjScoreGenes.internal <- function(
  codex_mRNA,
  gene_pairs,
  adj_matrix,
  c=0,
  num_cores=1,
  num_perms=1000,
  perm_estimate=T
) {
  codex_mRNA_cut <- t(codex_mRNA[,colnames(codex_mRNA) %in% as.vector(gene_pairs)])
  adjacency_score(adj_matrix, codex_mRNA_cut, gene_pairs, c=c, num_cores=num_cores,
                  num_perms=num_perms, perm_estimate = perm_estimate)
}
