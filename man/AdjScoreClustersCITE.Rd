% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_scores.R
\name{AdjScoreClustersCITE}
\alias{AdjScoreClustersCITE}
\title{Use the Adjacency Score to evaluate colocaliztion of all pairs
of CITE-seq clusters mapped to the CODEX spatial positions}
\usage{
AdjScoreClustersCITE(stvea_object, k, c = 0, num_cores = 1,
  num_perms = 1000, perm_estimate = T)
}
\arguments{
\item{stvea_object}{STvEA.data class object containing CITE-seq
cluster assignments, CODEX spatial coordinates, and transfer matrix
from GetTransferMatrix}

\item{k}{number of nearest neighbors to create graph
from CODEX spatial information}

\item{c}{constant used to determine width of diffusion, must be 0 <= c}

\item{num_cores}{integer specifying the number of cores to be used
in the computation. By default only one core is used.
On Windows, this must be set to 1.}

\item{num_perms}{number of permutations used to build the null
distribution for each feature. By default is set to 1000.}

\item{perm_estimate}{boolean indicating whether Gaussian distribution
parameters should be determined from num_perms permutations to estimate
the p-value. By default is set to TRUE.}
}
\description{
Use the Adjacency Score to evaluate colocaliztion of all pairs
of CITE-seq clusters mapped to the CODEX spatial positions
}
