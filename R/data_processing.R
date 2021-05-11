
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
#' @param model "nb" (Negative Binomial) or "gaussian" model to fit or
#' "arcsinh" for standard arcsinh transform (wrapper of flowCore method)
#' @param num_cores number of cores to use for parallelized fits
#' @param maxit maximum number of iterations for optim function
#'  - only used if model is "nb"
#' @param factr accuracy of optim function
#'  - only used if model is "nb"
#' @param optim_inits a matrix of (proteins x params) with initialization
#' parameters for each protein to input to the optim function. If NULL,
#' starts at two default parameter sets and picks the better one
#'  - only used if model is "nb"
#' @param normalize divide cleaned CODEX expression by total expression per cell
#'  - only used if model is "nb"
#'
#' @importFrom flowCore arcsinhTransform transformList transform
#'
#' @export
#'
CleanCODEX <- function(stvea_object,
                       model = "gaussian",
                       num_cores = 1,
                       maxit = 500,
                       factr = 1e-9,
                       optim_inits = NULL,
                       normalize = FALSE) {
  if (model == "gaussian") {
    # Normalize by total counts per cell
    codex_protein_norm <- stvea_object@codex_protein - min(stvea_object@codex_protein)
    avg_cell_total <- mean(rowSums(codex_protein_norm)) # multiply by constant to have nice numbers
    codex_protein_norm <- NormCells(codex_protein_norm) * avg_cell_total

    stvea_object@codex_clean <- CleanCODEX.gaussian.internal(codex_protein_norm)
  } else if (model == "nb") {
    stvea_object@codex_clean <- CleanCODEX.nb.internal(stvea_object@codex_protein,
                                                     num_cores=num_cores,
                                                     factr=factr,
                                                     optim_inits=optim_inits,
                                                     normalize=normalize)
  } else if (model == "arcsinh") {
    asinhTrans <- arcsinhTransform(a=0, b=1/5)
    translist <- transformList(colnames(stvea_object@codex_protein), asinhTrans)
    codex_protein_clean <- transform(stvea_object@codex_protein, translist)
    # Normalize each protein between [0,1]
    codex_protein_clean <- t(codex_protein_clean) - apply(codex_protein_clean,2,min)
    codex_protein_clean <- t(codex_protein_clean/ apply(codex_protein_clean,1,max))
    stvea_object@codex_clean <- codex_protein_clean
  } else {
    stop(sprintf("Invalid model parameter %s for CleanCODEX",model))
  }
  return(stvea_object)
}


#' Removes noise from CITE-seq protein data by fitting
#' a two component mixture model and computing each
#' expression measurement as the cumulative distribution
#' of the component with the higher median.
#'
#' This mixture model can either be:
#' - a Negative Binomial on the
#' expression counts (with optional weighting/normalization by
#' the total ADT counts per cell after cleaning)
#' - a Gaussian on the log-normalized expression with
#' zeros removed, similar to the method proposed by Trong et al. in SISUA
#' (https://www.biorxiv.org/content/10.1101/631382v1)
#'
#' @param stvea_object STvEA.data class with CITE-seq protein data
#' @param model "nb" (Negative Binomial) or "gaussian" model to fit
#' @param num_cores number of cores to use for parallelized fits
#' @param maxit maximum number of iterations for optim function
#'  - only used if model is "nb"
#' @param factr accuracy of optim function
#'  - only used if model is "nb"
#' @param optim_inits a matrix of (proteins x params) with initialization
#' parameters for each protein to input to the optim function. If NULL,
#' starts at two default parameter sets and picks the better one
#'  - only used if model is "nb"
#' @param normalize divide cleaned CITE-seq expression by total ADT counts per cell
#'  - only used if model is "nb"
#' @export
#'
CleanCITE <- function(stvea_object,
                      model = "nb",
                      num_cores = 1,
                      maxit = 500,
                      factr = 1e-9,
                      optim_inits = NULL,
                      normalize = TRUE) {
  if (is.null(stvea_object@cite_protein)) {
    stop("Input object must contain CITE-seq protein expression", call. =FALSE)
  }
  if (model == "nb") {
    stvea_object@cite_clean <- CleanCITE.nb.internal(stvea_object@cite_protein,
                                                  num_cores=num_cores,
                                                  factr=factr,
                                                  optim_inits=optim_inits,
                                                  normalize=normalize)
  } else if (model == "gaussian") {
    norm_protein <- NormCells(stvea_object@cite_protein)
    log_norm_protein <- log(1 + 1e4*norm_protein)
    stvea_object@cite_clean <- CleanCITE.gaussian.internal(log_norm_protein,
                                                           num_cores=num_cores)
  } else {
    stop(sprintf("Invalid model parameter %s for CleanCITE",model))
  }
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
  # Create filters
  if (is.null(size_lim)) {
    size_lim <- quantile(size, probs=c(0.025,0.99))
  }
  size_filter <- size >= size_lim[1] & size <= size_lim[2]

  if (is.null(blank_upper)) {
    blank_upper <- apply(blanks, 2, quantile, probs=0.995)
  }
  if (is.null(blank_lower)) {
    blank_lower <- apply(blanks, 2, quantile, probs=0.002)
  }
  blank_filter <- rep(FALSE, nrow(blanks))
  for (i in 1:ncol(blanks)) {
    blank_filter <- blank_filter | (blanks[,i] >= blank_lower[i])
  }
  for (i in 1:ncol(blanks)) {
    blank_filter <- blank_filter & (blanks[,i] <= blank_upper[i])
  }

  # Apply filters
  codex_filtered <- codex_raw[size_filter & blank_filter,]
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
CleanCODEX.gaussian.internal <- function(codex_filtered) {
  codex_clean <- codex_filtered
  # For each protein
  for (i in 1:ncol(codex_filtered)) {
    # Compute Gaussian mixture on each protein
    fit = Mclust(codex_filtered[,i], G=2, model="V", verbose=FALSE)
    signal <- as.numeric(which.max(fit$parameters$mean))
    # Compute cleaned data from cumulative of higher mean Gaussian
    expr_clean <- pnorm(codex_filtered[,i], mean = fit$parameters$mean[signal],
                        sd = sqrt(fit$parameters$variance$sigmasq[signal]))
    codex_clean[,i] <- expr_clean
  }
  return(codex_clean)
}

#' Removes noise from CODEX protein data by fitting
#' a two-component Negative Binomial mixture and computing
#' each expression measurement as the cumulative distribution
#' of the Negative Binomial with the higher median.
#' Normalizes CODEX protein data by the original
#' total counts per cell.
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param codex_protein Raw CODEX protein data (cell x protein)
#' @param num_cores number of cores to use for parallelized fits.
#' On Windows, this must be set to 1.
#' @param maxit maximum number of iterations for optim function
#' @param factr accuracy of optim function
#' @param optim_inits a matrix of (proteins x params) with initialization
#' parameters for each protein to input to the optim function. If NULL,
#' starts at two default parameter sets and picks the better one
#' @param normalize divide cleaned CODEX expression by total expression per cell
#'
#' @return Cleaned CODEX protein data matrix (cell x protein)
#'
#' @importFrom parallel mclapply
#'
#' @export
#'
CleanCODEX.nb.internal <- function(codex_protein,
                                  num_cores = 1,
                                  maxit = 500,
                                  factr = 1e-9,
                                  optim_inits = NULL,
                                  normalize = FALSE) {
  codex_protein_ceil <- ceiling(codex_protein) # NB requires integer data
  # could use round() for this, but CyTOF randomizes counts with value in [-1,0] so ceiling works for both

  if (!is.null(optim_inits)) {
    if (num_cores > 1) {
      codex_protein_list <- mclapply(1:ncol(codex_protein_ceil),
                                    function(i) FitNB(codex_protein_ceil[,i],
                                                      maxit=maxit, factr=factr,
                                                      optim_init = optim_inits[i,]),
                                    mc.cores=num_cores)
    } else {
      codex_protein_list <- lapply(1:ncol(codex_protein_ceil),
                                  function(i) FitNB(codex_protein_ceil[,i],
                                                    maxit=maxit, factr=factr,
                                                    optim_init = optim_inits[i,]))
    }
  } else {
    if (num_cores > 1) {
      codex_protein_list <- mclapply(1:ncol(codex_protein_ceil),
                                    function(i) FitNB(codex_protein_ceil[,i],
                                                      maxit=maxit, factr=factr),
                                    mc.cores=num_cores)
    } else {
      codex_protein_list <- lapply(1:ncol(codex_protein_ceil),
                                  function(i) FitNB(codex_protein_ceil[,i],
                                                    maxit=maxit, factr=factr))
    }
  }

  codex_protein_clean <- codex_protein
  for (i in 1:ncol(codex_protein_clean)) {
    codex_protein_clean[,i] <- codex_protein_list[[i]]
  }

  if (normalize) {
    # Normalize by original total expression per cell
    codex_protein_clean <- NormCells(codex_protein_clean, codex_protein)
    codex_protein_clean <- t(codex_protein_clean) - apply(codex_protein_clean,2,min)
    codex_protein_clean <- t(codex_protein_clean/ apply(codex_protein_clean,1,max))
  }
  return(codex_protein_clean)
}


#' Removes noise from CITE-seq protein data by fitting
#' a two-component Gaussian mixture to the log-normalized
#' expression with the zeros removed, and then computing each
#' expression measurement as the cumulative distribution
#' of the Gaussian with the higher median.
#'
#' This model is similar to the one proposed by Trong et al. in the SISUA preprint
#' (https://www.biorxiv.org/content/10.1101/631382v1)
#'
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param norm_cite_protein Log-normalized CITE-seq protein data (cell x protein)
#' @param num_cores number of cores to use for parallelized fits.
#' On Windows, this must be set to 1.
#'
#' @return Cleaned CITE-seq protein data matrix (cell x protein)
#'
#' @importFrom parallel mclapply
#'
#' @export
#'
CleanCITE.gaussian.internal <- function(norm_cite_protein, num_cores = 1) {
  if (num_cores > 1) {
    cite_protein_list <- mclapply(1:ncol(norm_cite_protein),
                                  function(i) FitGaussian(norm_cite_protein[,i]),
                                  mc.cores=num_cores)
  } else {
    cite_protein_list <- lapply(1:ncol(norm_cite_protein),
                                function(i) FitGaussian(norm_cite_protein[,i]))
  }

  cite_protein_clean <- norm_cite_protein
  for (i in 1:ncol(norm_cite_protein)) {
    cite_protein_clean[,i] <- cite_protein_list[[i]]
  }

  return(cite_protein_clean)
}


#' Removes noise from CITE-seq protein data by fitting
#' a two-component Negative Binomial mixture and computing
#' each expression measurement as the cumulative distribution
#' of the Negative Binomial with the higher median.
#' Normalizes CITE-seq protein data by the original
#' total counts per cell.
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param cite_protein Raw CITE-seq protein data (cell x protein)
#' @param num_cores number of cores to use for parallelized fits.
#' On Windows, this must be set to 1.
#' @param maxit maximum number of iterations for optim function
#' @param factr accuracy of optim function
#' @param optim_inits a matrix of (proteins x params) with initialization
#' parameters for each protein to input to the optim function. If NULL,
#' starts at two default parameter sets and picks the better one
#' @param normalize divide cleaned CITE-seq expression by total ADT counts per cell
#'
#' @return Cleaned CITE-seq protein data matrix (cell x protein)
#'
#' @importFrom parallel mclapply
#'
#' @export
#'
CleanCITE.nb.internal <- function(cite_protein,
                               num_cores = 1,
                               maxit = 500,
                               factr = 1e-9,
                               optim_inits = NULL,
                               normalize = TRUE) {
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

  if (normalize) {
    # Normalize by original total counts per cell
    cite_protein_clean <- NormCells(cite_protein_clean, cite_protein)
    cite_protein_clean <- t(cite_protein_clean) - apply(cite_protein_clean,2,min)
    cite_protein_clean <- t(cite_protein_clean/ apply(cite_protein_clean,1,max))
  }
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
    # Sometimes negative binomial doesn't fit well with certain starting parameters, so try 2
    fit1 <- optim(c(5,50,2,0.5,0.5), SSE, p_obs = p_obs,
                  method="L-BFGS-B", lower=rep(1e-8,5), upper=c(Inf,Inf,Inf,Inf,1),
                  control=list(maxit=maxit, factr=factr))
    fit2 <- optim(c(5,50,0.5,2,0.5), SSE, p_obs = p_obs,
                  method="L-BFGS-B", lower=rep(1e-8,5), upper=c(Inf,Inf,Inf,Inf,1),
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
                 method="L-BFGS-B", lower=rep(1e-8,5), upper=c(Inf,Inf,Inf,Inf,1),
                 control=list(maxit=maxit, factr=factr))$par
  }

  # Distribution with higher median is signal
  signal <- as.numeric(which.max(c(qnbinom(0.5,mu=fit[1],size=1/fit[3]),
                                   qnbinom(0.5,mu=fit[2],size=1/fit[4]))))

  expr_clean <- pnbinom(protein_expr, mu=fit[signal], size=1/fit[signal+2])
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


#' Fits the log-normalized expression values of one protein with a Gaussian mixture
#' Takes matrices and data frames instead of STvEA.data class
#'
#' @param norm_protein_expr Log-normalized CITE-seq protein data for one protein
#'
#' @importFrom mclust Mclust
#'
FitGaussian <- function(norm_protein_expr) {
  npe <- norm_protein_expr[norm_protein_expr != 0]
  fit = Mclust(npe, G=2, model="V", verbose=FALSE)
  signal <- as.numeric(which.max(fit$parameters$mean))
  expr_clean <- pnorm(norm_protein_expr,
                 mean = fit$parameters$mean[signal],
                 sd = sqrt(fit$parameters$variance$sigmasq[signal]))
  return(expr_clean)
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

#' Normalize protein matrix by dividing out total counts per cell,
#' accounting for cells with no counts by keeping them 0
#' Usually to_norm and norm_by will be the same count matrix,
#' but sometimes we want to normalize processed data by raw
#'
#' @param to_norm a (cells x proteins) matrix
#' @param norm_by a (cells x proteins) matrix
#'
NormCells <- function(to_norm, norm_by=to_norm) {
  nonzero <- rowSums(norm_by) != 0
  protein_norm <- to_norm
  protein_norm[nonzero,] <- protein_norm[nonzero,] / rowSums(norm_by[nonzero,])
  return(protein_norm)
}

