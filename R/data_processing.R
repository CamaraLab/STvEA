
# Function wrappers using STvEA.data class

#' Removes points from CODEX matrix that are not cells
#' as determined by the gating strategy on the blank channels
#' from the CODEX paper
#'
#' @param stvea_object STvEA.data class with CODEX data
#' @param size_lim lower and upper limits on size of each cell. If blank, set to 0.025 and 0.99 quantiles
#' @param blank_upper a vector with an upper bound expression cutoff for each blank channel.
#' If NULL, blank upper bounds are set as the 0.995 quantile for each blank
#' @param blank_lower a vector with a lower bound expression cutoff for each blank channel.
#' If NULL, blank lower bounds are set as the 0.002 quantile for each blank
#'
#' @export
#'
FilterCODEX <- function(stvea_object,
                                size_lim = NULL,
                                blank_upper = NULL,
                                blank_lower = NULL) {
  if (is.null(stvea_object@codex_protein) || is.null(stvea_object@codex_size) ||
      is.null(stvea_object@codex_blanks)) {
    stop("Input object must contain size of each CODEX cell and expression for CODEX protein and blank channels", call. =FALSE)
  }
  if (is.null(row.names(stvea_object@codex_protein))) {
    row.names(stvea_object@codex_protein) <- as.character(1:nrow(stvea_object@codex_protein))
  }
  filter_matrix <- FilterCODEX.internal(stvea_object@codex_protein,
                                                      stvea_object@codex_size,
                                                      stvea_object@codex_blanks,
                                                      size_lim=size_lim,
                                                      blank_upper=blank_upper,
                                                      blank_lower=blank_lower)
  filter <- row.names(stvea_object@codex_protein) %in% row.names(filter_matrix)
  stvea_object@codex_protein <- filter_matrix
  stvea_object@codex_size <- stvea_object@codex_size[filter]
  stvea_object@codex_blanks <- stvea_object@codex_blanks[filter,]

  if (!is.null(stvea_object@codex_spatial)) {
    stvea_object@codex_spatial <- stvea_object@codex_spatial[filter,]
  }
  if (!is.null(stvea_object@codex_clusters)) {
    stvea_object@codex_clusters <- stvea_object@codex_clusters[filter]
  }
  if (!is.null(stvea_object@codex_emb)) {
    stvea_object@codex_emb <- stvea_object@codex_emb[filter,]
  }
  return(stvea_object)
}


#' Removes noise from CODEX data by fitting a Gaussian
#' mixture and computing each expression measurement
#' as the cumulative distribution of the Gaussian with
#' the higher mean.
#'
#' @param stvea_object STvEA.data class with CODEX data after FilterCODEX
#'
#' @export
#'
CleanCODEX <- function(stvea_object) {
  stvea_object@codex_clean <- CleanCODEX.internal(stvea_object@codex_protein)
  return(stvea_object)
}


#' Removes noise from CITE-seq protein data by fitting
#' a Negative Binomial mixture and computing each expression
#' measurement as the cumulative distribution of the
#' Negative Binomial with the higher median.
#'
#' @param stvea_object STvEA.data class with CITE-seq protein data
#' @param maxit maximum number of iterations for optim function
#' @param factr accuracty of optim function
#' @param optim_inits a matrix of (proteins x params) with initialization
#' parameters for each protein to input to the optim function. If NULL,
#' starts at two default parameter sets and picks the better one
#' @param num_cores number of cores to use for parallelized fits
#'
#' @export
#'
CleanCITE <- function(stvea_object,
                      maxit = 500,
                      factr = 1e-9,
                      optim_inits = NULL,
                      num_cores = 1) {
  if (is.null(stvea_object@cite_protein)) {
    stop("Input object must contain CITE-seq protein expression", call. =FALSE)
  }
  stvea_object@cite_clean <- CleanCITE.internal(stvea_object@cite_protein,
                                                  factr=factr,
                                                  optim_inits=optim_inits,
                                                  num_cores=num_cores)
  return(stvea_object)
}


#' Filter CITE-seq mRNA data based on number of transcripts per cell
#' and the number of cells each gene is expressed in
#'
#' @param stvea_object STvEA.data class with CITE-seq mRNA data
#' @param min_transcripts minimum number of transcripts to keep cell
#' @param min_cells_per_gene minimum number of cells a gene must be
#'
#' @export
#'
FilterCITEmRNA <- function(stvea_object, min_transcripts=1200, min_cells_per_gene=30) {
  if (is.null(stvea_object@cite_mRNA)) {
    stop("stvea_object must contain CITE-seq mRNA data", call. =FALSE)
  }
  stvea_object@cite_mRNA <- FilterCITEmRNA.internal(stvea_object@cite_mRNA,
                                                    min_transcripts = min_transcripts,
                                                    min_cells_per_gene = min_cells_per_gene)
}


# Functions with matrix parameters, not using STvEA.data class

#' Removes points from CODEX matrix that are not cells
#' as determined by the gating strategy on the blank channels
#' from the CODEX paper, then normalizes data by total counts per cell
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param codex_raw CODEX expression matrix after spillover correction,
#' not including blank channels (cells x proteins)
#' @param size vector of cell sizes from CODEX segmentation
#' @param blanks expression matrix for blank channels after spillover correction (cells x channels)
#' @param size_lim lower and upper limits on size of each cell.
#' If blank, set to 0.025 and 0.99 quantiles
#' @param blank_upper a vector with an upper bound expression cutoff for each blank channel.
#' If NULL, blank upper bounds are set as the 0.995 quantile for each blank
#' @param blank_lower a vector with a lower bound expression cutoff for each blank channel.
#' If NULL, blank lower bounds are set as the 0.002 quantile for each blank
#'
#' @return CODEX expression matrix with some cells filtered out (cell x protein)
#'
#' @export
#'
FilterCODEX.internal <- function(codex_raw, size, blanks,
                         size_lim = NULL,
                         blank_upper = NULL,
                         blank_lower = NULL) {
  # filter cells by size
  # determine cutoffs based on blank expression, but also let user define
  # plot blank channels to user to verify
  # normalize by total expression per cell
  # return filtered matrix (maybe as part of object)

  # Create filters
  if (is.null(size_lim)) {
    size_lim <- quantile(size, probs=c(0.025,0.99))
  }
  size_filter <- size > size_lim[1] & size < size_lim[2]

  if (is.null(blank_upper)) {
    blank_upper <- apply(blanks, 2, quantile, probs=0.995)
  }
  if (is.null(blank_lower)) {
    blank_lower <- apply(blanks, 2, quantile, probs=0.002)
  }
  blank_filter <- rep(FALSE, nrow(blanks))
  for (i in 1:ncol(blanks)) {
    blank_filter <- blank_filter | (blanks[,i] > blank_lower[i])
  }
  for (i in 1:ncol(blanks)) {
    blank_filter <- blank_filter & (blanks[,i] < blank_upper[i])
  }

  # Normalize by total counts per cell
  codex_filtered <- codex_raw - min(codex_raw)
  avg_cell_total <- mean(rowSums(codex_filtered)) # multiply by constant to have nice numbers
  codex_filtered <- codex_filtered/rowSums(codex_filtered) * avg_cell_total

  # Apply filters
  codex_filtered <- codex_filtered[size_filter & blank_filter,]
  return(codex_filtered)
}


#' Removes noise from CODEX data by fitting a Gaussian
#' mixture and computing each expression measurement
#' as the cumulative distribution of the Gaussian with
#' the higher mean.
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param codex_filtered CODEX expression matrix after filtering (cell x proteins)
#'
#' @return Cleaned CODEX protein expression matrix (cell x protein)
#'
#' @import mclust
#'
#' @export
#'
CleanCODEX.internal <- function(codex_filtered) {
  # compute Gaussian mixture on each protein
  # compute cleaned data from cumulative of higher mean Gaussian
  # return cleaned data (maybe as part of object)

  codex_clean <- codex_filtered
  for (i in 1:ncol(codex_filtered)) {
    fit = Mclust(codex_filtered[,i], G=2, model="V", verbose=FALSE)
    signal <- as.numeric(which.max(fit$parameters$mean))
    expr_clean <- pnorm(codex_filtered[,i], mean = fit$parameters$mean[signal],
                        sd = sqrt(fit$parameters$variance$sigmasq[signal]))
    codex_clean[,i] <- expr_clean
  }
  return(codex_clean)
}


#' Removes noise from CITE-seq protein data by fitting
#' a Negative Binomial mixture and computing each expression
#' measurement as the cumulative distribution of the
#' Negative Binomial with the higher median.
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param cite_protein Raw CITE-seq protein data (cell x protein)
#' @param maxit maximum number of iterations for optim function
#' @param factr accuracty of optim function
#' @param optim_inits a matrix of (proteins x params) with initialization
#' parameters for each protein to input to the optim function. If NULL,
#' starts at two default parameter sets and picks the better one
#' @param num_cores number of cores to use for parallelized fits.
#' On Windows, this must be set to 1.
#'
#' @return Cleaned CITE-seq protein data matrix (cell x protein)
#'
#' @importFrom parallel mclapply
#'
#' @export
#'
CleanCITE.internal <- function(cite_protein,
                       maxit = 500,
                       factr = 1e-9,
                       optim_inits = NULL,
                       num_cores = 1) {
  # Fit Negative Binomial mixture to protein data
  # Calculate cleaned data from cumulative of higher median

  if (!is.null(optim_inits)) {
    if (num_cores > 1) {
      cite_protein_list <- mclapply(1:ncol(cite_protein),
                                    function(i) FitNB(cite_protein[,i],
                                                            maxit=maxit, factr=factr,
                                                            optim_init = optim_inits[i,]),
                                    mc.cores=num_cores)
    } else {
      cite_protein_list <- lapply(1:ncol(cite_protein),
                                    function(i) FitNB(cite_protein[,i],
                                                      maxit=maxit, factr=factr,
                                                      optim_init = optim_inits[i,]))
    }
  } else {
    if (num_cores > 1) {
      cite_protein_list <- mclapply(1:ncol(cite_protein),
                                    function(i) FitNB(cite_protein[,i],
                                                            maxit=maxit, factr=factr),
                                    mc.cores=num_cores)
    } else {
      cite_protein_list <- lapply(1:ncol(cite_protein),
                                    function(i) FitNB(cite_protein[,i],
                                                      maxit=maxit, factr=factr))
    }
  }

  cite_protein_clean <- cite_protein
  for (i in 1:ncol(cite_protein)) {
    cite_protein_clean[,i] <- cite_protein_list[[i]]
  }

  cite_protein_clean <- cite_protein_clean / rowSums(cite_protein)
  cite_protein_clean <- t(cite_protein_clean) - apply(cite_protein_clean,2,min)
  cite_protein_clean <- t(cite_protein_clean/ apply(cite_protein_clean,1,max))
  return(cite_protein_clean)
}


#' Fits the expression values of one protein with a Negative Binomial mixture
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param protein_expr Raw CITE-seq protein data for one protein
#' @param maxit maximum number of iterations for optim function
#' @param factr accuracty of optim function
#' @param optim_init optional initialization parameters for the optim function
#' if NULL, starts at two default parameter sets and picks the better one
#'
#' @importFrom stats optim
#'
FitNB <- function(protein_expr,
                        maxit = 500,
                        factr = 1e-9,
                        optim_init = NULL) {
  p_obs <- table(factor(protein_expr, levels=0:max(protein_expr)))/length(protein_expr)

  if (is.null(optim_init)) {
    # sometimes negative binomial doesn't fit well with certain starting parameters, so try 2
    fit1 <- optim(c(5,50,2,0.5,0.5), SSE, p_obs = p_obs,
                  method="L-BFGS-B", lower=rep(1e-8,5), upper=c(Inf, Inf,Inf,Inf,1),
                  control=list(maxit=maxit, factr=factr))
    fit2 <- optim(c(5,50,0.5,2,0.5), SSE, p_obs = p_obs,
                  method="L-BFGS-B", lower=rep(1e-8,5), upper=c(Inf, Inf,Inf,Inf,1),
                  control=list(maxit=maxit, factr=factr))
    score1 <- SSE(fit1$par, p_obs)
    score2 <- SSE(fit2$par, p_obs)
    if (score1 < score2) {
      fit <- fit1$par
    } else {
      fit <- fit2$par
    }
  } else {
    fit <- optim(optim_init, SSE, p_obs = p_obs,
                 method="L-BFGS-B", lower=rep(1e-8,5), upper=c(Inf, Inf,Inf,Inf,1),
                 control=list(maxit=maxit, factr=factr))$par
  }

  # distribution with higher median is signal
  signal <- as.numeric(which.max(c(qnbinom(0.5,mu=fit[1],size=1/fit[3]),
                                   qnbinom(0.5,mu=fit[2],size=1/fit[3]))))

  expr_clean <- pnbinom(protein_expr, mu=fit[signal],size=1/fit[signal+2])
  return(expr_clean)
}


#' Calculates the sum of squared errors in binned probabilities of count data
#'
#' @param args arguments used in the negative binomial mixture model
#' @param p_obs a named vector of the probabilities of observing a given count
#'  in gene expression data, as output by running table() on the gene count data
#'
SSE <- function(args, p_obs) {
  p_exp <- args[5]*dnbinom(as.numeric(names(p_obs)), mu=args[1], size=1/args[3]) +
    (1-args[5])*dnbinom(as.numeric(names(p_obs)), mu=args[2], size=1/args[4])
  sse <- min(sum((p_exp - p_obs)^2), .Machine$integer.max)
  return(sse)
}


#' Filter CITE-seq mRNA data based on number of transcripts per cell
#' and the number of cells each gene is expressed in
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param cite_mRNA a (cell x gene) matrix of CITE-seq mRNA data
#' @param min_transcripts minimum number of transcripts to keep cell
#' @param min_cells_per_gene minimum number of cells a gene must be
#' expressed in to be kept
#'
FilterCITEmRNA.internal <- function(cite_mRNA,
                                    min_transcripts = 1200,
                                    min_cells_per_gene = 30) {
  transcripts <- rowSums(cite_mRNA)
  cite_mRNA <- cite_mRNA[transcripts >= min_transcripts,]
  cells_per_gene <- colSums(f != 0)
  cite_mRNA <- cite_mRNA[cells_per_gene >= min_cells_per_gene,]
  return(cite_mRNA)
}


