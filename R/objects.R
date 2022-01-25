
#' Class to hold all data matrices needed for mapping CODEX to CITE-seq
#' and following analysis
#'
STvEA.data <- setClass(
  Class = 'STvEA.data',
  slots = c(
    cite_mRNA = 'ANY', # cite_cells x genes
    cite_mRNA_norm = 'ANY', # cite cells x genes - log(1+TPM) or seurat
    cite_latent = 'ANY', # cite cells x low
    cite_emb = 'ANY', # cite cells x 2
    hdbscan_param_scan = 'list', # output of ParameterScan
    cite_clusters = 'vector',
    cite_protein = 'ANY', # cite cells x proteins
    cite_clean = 'ANY', # cite cells x proteins
    cite_norm = 'ANY', # cite cells x proteins

    codex_protein = 'ANY', # codex cells x proteins
    codex_size = 'numeric',
    codex_blanks = 'ANY', # codex cells x blank channels
    codex_spatial = 'ANY', # codex cells x 3
    codex_emb = 'ANY', # codex cells x 2
    codex_knn = 'ANY', # codex cells x k
    codex_clusters = 'vector',
    codex_clean = 'ANY', # codex cells x proteins

    corrected_codex = 'ANY', # codex cells x proteins
    transfer_matrix = 'dgCMatrix', # codex cells x cite cells
                                 # k CODEX neighbors for each CITE-seq cell
    codex_mRNA = 'ANY' # codex cells x genes
  )
)


#' Set CODEX data in STvEA.data object
#'
#' @param codex_protein Segmented and spillover corrected CODEX protein expression (cell x protein)
#' @param codex_blanks Segmented and spillover corrected CODEX blank channel expression (cell x blank channel)
#' @param codex_size Size of each CODEX cell (vector)
#' @param codex_spatial (optional) xyz coordinates of each CODEX cell (cell x 3 coordinates)
#' @param stvea_object (optional) Pre-existing STvEA.data object to load data into.
#' If not provided, a new object is created.
#'
#' @export
#'
SetDataCODEX <- function(codex_protein,
                         codex_blanks,
                         codex_size,
                         codex_spatial = NULL,
                         stvea_object = NULL) {
  if (is.null(colnames(codex_protein))) {
    warning("No protein names in CODEX protein matrix. Assigning names CODEXprotein1, CODEXprotein2, ...")
    colnames(codex_protein) <- paste("CODEXprotein",1:ncol(codex_protein),sep="")
  }
  if (is.null(stvea_object)) {
    stvea_object <- new(
      Class = "STvEA.data",
      codex_protein = codex_protein,
      codex_blanks = codex_blanks,
      codex_size = codex_size,
      codex_spatial = codex_spatial
    )
  } else {
    stvea_object@codex_protein <- codex_protein
    stvea_object@codex_blanks <- codex_blanks
    stvea_object@codex_size <- codex_size
    stvea_object@codex_spatial <- codex_spatial
  }
  return(stvea_object)
}


#' Set CITE-seq data in STvEA.data object
#'
#' @param cite_mRNA Raw count data for CITE-seq mRNA (cell x gene)
#' @param cite_protein Raw expression data for CITE-seq proteins (cell x protein)
#' @param cite_mRNA_norm (optional) Normalized and scaled CITE-seq mRNA count data. (cell x gene)
#' If not provided, raw count data is divided by total per cell, then log transformed.
#' @param cite_latent (optional) Low dimensional latent space for CITE-seq mRNA (cell x dim)
#' @param stvea_object (optional) Pre-existing STvEA.data object to load data into.
#' If not provided, a new object is created.
#'
#' @export
#'
SetDataCITE <- function(cite_mRNA,
                        cite_protein,
                        cite_mRNA_norm = NULL,
                        cite_latent = NULL,
                        stvea_object = NULL) {
  if (is.null(colnames(cite_protein))) {
    warning("No protein names in CITE-seq protein matrix. Assigning names CITEprotein1, CITEprotein2, ...")
    colnames(cite_protein) <- paste("CITEprotein",1:ncol(cite_protein),sep="")
  }
  if (is.null(cite_mRNA_norm)) {
    cite_mRNA_norm <- log(1 + 1e4*(cite_mRNA/rowSums(cite_mRNA)))
  }
  if (is.null(stvea_object)) {
    stvea_object <- new(
      Class = "STvEA.data",
      cite_mRNA = cite_mRNA,
      cite_mRNA_norm = cite_mRNA_norm,
      cite_protein = cite_protein
    )
  } else {
    stvea_object@cite_mRNA <- cite_mRNA
    stvea_object@cite_mRNA_norm <- cite_mRNA_norm
    stvea_object@cite_protein <- cite_protein
  }
  if(!is.null(cite_latent)) {
    stvea_object@cite_latent <- cite_latent
  }
  return(stvea_object)
}


#' Transfer data from a Seurat object to STvEA.data object
#'
#' @param seurat_object Seurat class object
#' @param embedding_reduction (optional) name of the slot holding the dimension reduction
#' embedding used for visualization. Examples: "umap", "tsne"
#' @param latent_reduction (optional) name of the lot holding the dimension reduction to
#' be used for the CITE-seq mRNA latent space. Examples: "pca"
#' @param latent_dims number of dimensions to use for the latent space,
#' not used if latent_reduction is not defined
#'
#' @return STvEA.data class object
#'
#' @export
#'
TransferDataSeurat2 <- function(seurat_object,
                                embedding_reduction = NULL,
                                latent_reduction = NULL,
                                latent_dims = 20) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  stvea_object <- new(
    Class = "STvEA.data",
    cite_mRNA = t(seurat_object@raw.data),
    cite_protein = t(Seurat::GetAssayData(seurat_object, assay.type="CITE", slot="raw.data")),
    cite_clusters = Seurat::GetClusters(seurat_object)$cluster
  )
  if (!is.null(seurat_object@scale.data)) {
    stvea_object@cite_mRNA_norm <- t(seurat_object@scale.data)
  }
  if (!is.null(latent_reduction)) {
    cite_latent <- Seurat::GetDimReduction(
      object = seurat_object,
      reduction.type = latent_reduction,
      slot = 'cell.embeddings'
    )
    stvea_object@cite_latent <- cite_latent[,1:min(latent_dims,ncol(cite_latent))]
  }
  if (!is.null(embedding_reduction)) {
    stvea_object@cite_emb <- as.data.frame(Seurat::GetDimReduction(
      object = seurat_object,
      reduction.type = embedding_reduction,
      slot = 'cell.embeddings'
    ))
  }
  return(stvea_object)
}


#' Read in FCS file
#'
#' @param path_to_fcs readable path to FCS file
#' @param is_protein boolean vector indicating which columns in FCS
#' should be kept as protein expression for downstream analyses
#' @param is_blank boolean vector indicating which columns in FCS
#' are blank for filtering purposes
#' @param protein_names vector of names for each protein channel column -
#' if null, will keep column names from FCS
#' @param tile_relative whether x,y,z coordinates are relative to x_tile,y_tile (TRUE)
#' or to corner of image (FALSE).
#'       If TRUE, will compute absolute coordinates
#'       If FALSE, will just take given x,y,z coordinates
#' @param num_tile_x number of tiles in the x direction, only used if tile_relative = TRUE
#' and the tile information in the FCS file is saved as tile_nr instead of x_tile,y_tile
#' @param x_pix_size x dimension of pixel coordinates in nm. Set to 1 to keep pixel dimensions.
#' @param y_pix_size y dimension of pixel coordinates in nm. Set to 1 to keep pixel dimensions.
#' @param z_pix_size z dimension of pixel coordinates in nm. Set to 1 to keep pixel dimensions.
#' @param stvea_object (optional) Pre-existing STvEA.data object to load data into.
#' If not provided, a new object is created.
#'
#' @return STvEA.data class object
#'
#' @importFrom flowCore read.FCS
#' @importFrom stringr str_match
#'
#' @export
#'
ReadDataFCS <- function(path_to_fcs,
                        is_protein,
                        is_blank,
                        protein_names = NULL,
                        tile_relative = FALSE,
                        num_tile_x = 9,
                        x_pix_size = 188,
                        y_pix_size = 188,
                        z_pix_size = 900,
                        stvea_object = NULL) {
  expr_mat <- read.FCS(path_to_fcs,transformation=FALSE, truncate_max_range=FALSE)@exprs

  # Get blank columns
  if (sum(is_blank) == 0) {
    warning("No columns are labelled as blank, setting blank channel data to all 0s")
    codex_blanks <- as.matrix(rep(0,nrow(expr_mat)), ncol=1)
  } else {
    codex_blanks <- expr_mat[,is_blank,drop=FALSE]
  }

  # Get protein matrix
  if (!is.null(protein_names) && sum(is_protein) != length(protein_names)) {
    stop("Length of protein_names must be the same as TRUE values in is_protein")
  }
  codex_protein <- expr_mat[,is_protein,drop=FALSE]
  if (!is.null(protein_names)) {
    colnames(codex_protein) <- protein_names
  }

  # Match columns of the format "X.X" and convert them to just "x"
  reg_match <- str_match(colnames(expr_mat), "(.+)\\.\\1")
  match <- !is.na(reg_match[,1]) # these columns have the form "X.X"
  colnames(expr_mat)[match] <- reg_match[match,2] # convert to just "X"
  colnames(expr_mat) <- sapply(colnames(expr_mat), tolower) # lowercase "x"

  # Get spatial information
  if (!all(c("x","y","z") %in% colnames(expr_mat))) {
    warning("Cannot find x,y,z coordinates in FCS. Continuing without spatial information")
    codex_spatial_nm <- NULL
  } else {
    if (tile_relative) {
      if ("tile_nr" %in% colnames(expr_mat)) {
        x <- floor((expr_mat[,"tile_nr"]-1)/num_tile_x) * max(expr_mat[,"x"]) + expr_mat[,"x"]
        y <- ((expr_mat[,"tile_nr"] - 1) %% num_tile_x) * max(expr_mat[,"y"]) + expr_mat[,"y"]
      } else if ("tile_num" %in% colnames(expr_mat)) {
        x <- floor((expr_mat[,"tile_num"]-1)/num_tile_x) * max(expr_mat[,"x"]) + expr_mat[,"x"]
        y <- ((expr_mat[,"tile_num"] - 1) %% num_tile_x) * max(expr_mat[,"y"]) + expr_mat[,"y"]
      } else if ("x_tile" %in% colnames(expr_mat) && "y_tile" %in% colnames(expr_mat)) {
        x <- (expr_mat[,"x_tile"] - 1) * max(expr_mat[,"x"]) + expr_mat[,"x"]
        y <- (expr_mat[,"y_tile"] - 1) * max(expr_mat[,"y"]) + expr_mat[,"y"]
      } else {
        stop("Cannot find tile information in FCS to compute absolute x,y,z from relative")
      }
      z <- expr_mat[,"z"]
      codex_spatial <- cbind(x,y,z)
    } else {
      codex_spatial <- expr_mat[,c("x","y","z")]
    }
    codex_spatial_nm <- cbind(x=codex_spatial[,"x"]*x_pix_size,
                              y=codex_spatial[,"y"]*y_pix_size,
                              z=codex_spatial[,"z"]*z_pix_size)
  }

  # Get size information
  if ("size" %in% colnames(expr_mat)) {
    codex_size <- expr_mat[,"size"]
  } else {
    warning("Cannot find size information in FCS, setting all cell size to 0")
    codex_size <- rep(0, nrow(expr_mat))
  }

  stvea_object <- SetDataCODEX(codex_protein = codex_protein,
                              codex_blanks = codex_blanks,
                              codex_size = codex_size,
                              codex_spatial = codex_spatial_nm,
                              stvea_object = stvea_object)

  return(stvea_object)
}

#' Reads channel names from FCS two different ways
#'
#' @param path_to_fcs readable path to FCS file
#'
#' @return list where "columns" are the column names from the FCS data matrix,
#' and "channels" are the readable names in the description
#'
#' @importFrom flowCore read.FCS
#'
#' @export
#'
ReadNamesFCS <- function(path_to_fcs) {
  data <- read.FCS(path_to_fcs,transformation=FALSE, truncate_max_range=FALSE)
  columns <- colnames(data@exprs)
  names(columns) <- NULL
  channels <- as.character(data@description[paste("$P",1:ncol(data@exprs),"S",sep="")])
  return(list(columns = columns, channels=channels))
}
