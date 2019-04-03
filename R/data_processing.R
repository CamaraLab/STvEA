

#' Removes points from CODEX matrix that are not cells
#' as determined by the gating strategy on the blank channels
#' from the CODEX paper
#'
#' @param codex_raw CODEX expression matrix after spillover correction, not including blank channels (cells x proteins)
#' @param size vector of cell sizes from CODEX segmentation
#' @param blanks expression matrix for blank channels after spillover correction (cells x channels)
#' @param size_lim lower and upper limits on size of each cell. If blank, set to 0.025 and 0.99 quantiles
#' @param blank_upper a vector with an upper bound expression cutoff for each blank channel.
#' If NULL, blank upper bounds are set as the 0.995 quantile for each blank
#' @param blank_lower a vector with a lower bound expression cutoff for each blank channel.
#' If NULL, blank lower bounds are set as the 0.002 quantile for each blank
#'
#' @return CODEX expression matrix with some cells filtered out (cell x protein)
#'
#' @export
#'
filter_codex <- function(codex_raw, size, blanks,
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
#'
#' @param codex_filtered CODEX expression matrix after filtering (cell x proteins)
#'
#' @return Cleaned CODEX protein expression matrix (cell x protein)
#'
#' @import mclust
#'
#' @export
#'
clean_codex <- function(codex_filtered) {
  # compute Gaussian mixture on each protein
  # compute cleaned data from cumulative of higher mean Gaussian
  # return cleaned data (maybe as part of object)

  codex_clean <- codex_filtered
  for (i in 1:ncol(codex_filtered)) {
    fit = Mclust(codex_filtered[,i], G=2, model="V", verbose=FALSE)
    signal <- as.numeric(which.max(fit$parameters$mean))
    expr_clean <- pnorm(codex_filtered[,i], mean = fit$parameters$mean[signal], sd = sqrt(fit$parameters$variance$sigmasq[signal]))
    codex_clean[,i] <- expr_clean
  }
  return(codex_clean)
}


#' Removes noise from CITE-seq protein data by fitting
#' a Negative Binomial mixture and computing each expression
#' measurement as the cumulative distribution of the
#' Negative Binomial with the higher median.
#'
#' @param cite_protein Raw CITE-seq protein data (cell x protein)
#' @param maxit maximum number of iterations for optim function
#' @param factr accuracty of optim function
#' @param optim_inits a matrix of (proteins x params) with initialization
#' parameters for each protein to input to the optim function. If NULL,
#' starts at two default parameter sets and picks the better one
#' @param verbose print out each iteration of negative binomial fit loop
#'
#' @return Cleaned CITE-seq protein data matrix (cell x protein)
#'
#' @importFrom parallel mclapply
#'
#' @export
#'
clean_cite <- function(cite_protein,
                       maxit = 500,
                       factr = 1e-9,
                       optim_inits = NULL,
                       num_cores = 1) {
  # Fit Negative Binomial mixture to protein data
  # Calculate cleaned data from cumulative of higher median

  if (!is.null(optim_inits)) {
    cite_protein_list <- mclapply(1:ncol(cite_protein),
                                  function(i) fit_protein(cite_protein[,i],
                                                          maxit=maxit, factr=factr,
                                                          optim_init = optim_inits[i,]),
                                  mc.cores=num_cores)
  } else {
    cite_protein_list <- mclapply(1:ncol(cite_protein),
                                  function(i) fit_protein(cite_protein[,i],
                                                          maxit=maxit, factr=factr),
                                  mc.cores=num_cores)
  }
  cite_protein_clean <- do.call(cbind,cite_protein_list)
  return(cite_protein_clean)
}


#' Fits the expression values of one protein with a
#' Negative Binomial mixture
#'
#' @param protein_expr Raw CITE-seq protein data for one protein
#' @param maxit maximum number of iterations for optim function
#' @param factr accuracty of optim function
#' @param optim_init optional initialization parameters for the optim function
#' if NULL, starts at two default parameter sets and picks the better one
#'
#' @importFrom stats optim
#'
fit_protein <- function(protein_expr,
                        maxit = 500,
                        factr = 1e-9,
                        optim_init = NULL) {
  p_obs <- table(factor(protein_expr, levels=0:max(protein_expr)))/length(protein_expr)

  if (is.null(optim_init)) {
    # sometimes negative binomial doesn't fit well with certain starting parameters, so try 2
    fit1 <- optim(c(5,50,2,0.5,0.5), sse_fn, p_obs = p_obs,
                  method="L-BFGS-B", lower=rep(1e-8,5), upper=c(Inf, Inf,Inf,Inf,1),
                  control=list(maxit=maxit, factr=factr))
    fit2 <- optim(c(5,50,0.5,2,0.5), sse_fn, p_obs = p_obs,
                  method="L-BFGS-B", lower=rep(1e-8,5), upper=c(Inf, Inf,Inf,Inf,1),
                  control=list(maxit=maxit, factr=factr))
    score1 <- sse_fn(fit1$par, p_obs)
    score2 <- sse_fn(fit2$par, p_obs)
    if (score1 < score2) {
      fit <- fit1$par
    } else {
      fit <- fit2$par
    }
  } else {
    fit <- optim(optim_init, sse_fn, p_obs = p_obs,
                 method="L-BFGS-B", lower=rep(1e-8,5), upper=c(Inf, Inf,Inf,Inf,1),
                 control=list(maxit=maxit, factr=factr))$par
  }

  # distribution with higher median is signal
  signal <- as.numeric(which.max(c(qnbinom(0.5,mu=fit[1],size=1/fit[3]), qnbinom(0.5,mu=fit[2],size=1/fit[3]))))

  expr_clean <- pnbinom(protein_expr, mu=fit[signal],size=1/fit[signal+2])
  return(expr_clean)
}


#' Calculates the sum of squared errors in binned probabilities of count data
#'
#' @param args arguments used in the negative binomial mixture model
#' @param p_obs a named vector of the probabilities of observing a given count
#'  in gene expression data, as output by running table() on the gene count data
#'
sse_fn <- function(args, p_obs) {
  p_exp <- args[5]*dnbinom(as.numeric(names(p_obs)), mu=args[1], size=1/args[3]) +
    (1-args[5])*dnbinom(as.numeric(names(p_obs)), mu=args[2], size=1/args[4])
  sse <- min(sum((p_exp - p_obs)^2), .Machine$integer.max)
  return(sse)
}


#' Normalize CITE-seq ADT expression by original cell expression totals before cleaning
#' Then perform batch correction using the MNN method from Haghverdi et al.
#'
#' @param cite_protein_raw Raw CITE-seq protein data (cell x protein)
#' @param cite_protein_clean Cleaned CITE-seq protein data (cell x protein) output by clean_cite()
#' @param batch Batch indices for running MNN batch correction
#' @param k.batch k to use for MNN batch correction
#'
#' @return Normalized and batch corrected CITE-seq protein data matrix (cell x protein)
#'
#' @export
#'
normalize_cite <- function(cite_protein_clean,
                           cite_protein_raw,
                           batch = rep(1, nrow(cite_protein_raw)),
                           k.batch = 20) {
  # Normalize data by dividing out cell totals (if param TRUE)
  # Batch correction (if different batches)
  # Renormalize between 0 and 1?
  # Return cleaned data matrix (maybe as part of object?)

  # Normalize by total expression per cell
  cite_protein_clean <- cite_protein_clean / rowSums(cite_protein_raw)

  # Batch correction
  # batch_ids <- unique(batch)
  # if (length(batch_ids) > 1) {
  #   batch_data <- lapply(batch_ids, function(x) t(cite_protein_clean[batch==x,]))
  #
  #   corrected_data <- do.call(mnnCorrect, c(batch_data, k=k.batch))
  #
  #   cite_protein_corrected <- cite_protein_clean
  #   for (n in 1:length(batch_ids)) {
  #     batch_corrected <- t(corrected_data$corrected[[n]])
  #     colnames(batch_corrected) <- row.names(batch_data[[n]])
  #     row.names(batch_corrected) <- colnames(batch_data[[n]])
  #     cite_protein_corrected[batch == batch_ids[n],] <- batch_corrected
  #   }
  #   cite_protein_clean <- cite_protein_corrected
  # }
  cite_protein_clean <- t(cite_protein_clean) - apply(cite_protein_clean,2,min)
  cite_protein_clean <- t(cite_protein_clean/ apply(cite_protein_clean,1,max))
  return(cite_protein_clean)
}


#' Perform HDBSCAN consensus clustering on CITE-seq latent space
#'
cite_clustering <- function(cite_latent) {
}


#' Assign user input labels to CITE-seq clusters for future plotting
#'
clustering_names <- function(labels) {
}
