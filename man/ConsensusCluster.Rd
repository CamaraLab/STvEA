% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{ConsensusCluster}
\alias{ConsensusCluster}
\title{Keep only HDBSCAN results that pass a certain silhouette score cutoff
and create a dissimilarity matrix between cells from the clusterings
Perform agglomerative hierarchical clustering on the
consensus dissimilarity matrix}
\usage{
ConsensusCluster(stvea_object, silhouette_cutoff, inconsistent_value,
  min_cluster_size)
}
\arguments{
\item{stvea_object}{STvEA.data class object after ParameterScan has been run}

\item{silhouette_cutoff}{minimum silhouette score to keep clustering}

\item{inconsistent_value}{input parameter to python fcluster determining
where clusters are cut in the hierarchical tree}

\item{min_cluster_size}{cells in clusters smaller than this value are
assigned a cluster ID of -1, indicating no cluster assignment}
}
\description{
Keep only HDBSCAN results that pass a certain silhouette score cutoff
and create a dissimilarity matrix between cells from the clusterings
Perform agglomerative hierarchical clustering on the
consensus dissimilarity matrix
}
