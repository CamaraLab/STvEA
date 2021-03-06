% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_scores.R
\name{AdjScoreProteins}
\alias{AdjScoreProteins}
\title{Use the Adjacency Score to evaluate colocalization of given
pairs of proteins in the CODEX spatial dimensions}
\usage{
AdjScoreProteins(stvea_object, protein_pairs = NULL, k = 5, c = 0,
  num_cores = 1, num_perms = 1000, perm_estimate = TRUE)
}
\arguments{
\item{stvea_object}{STvEA.data class object containing CODEX
protein expression and CODEX spatial coordinates}

\item{protein_pairs}{a 2 column matrix of protein pairs where each row
specifies the names of the proteins in a pair. If NULL, all pairs
of columns in the codex_protein matrix are used.}

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
Use the Adjacency Score to evaluate colocalization of given
pairs of proteins in the CODEX spatial dimensions
}
