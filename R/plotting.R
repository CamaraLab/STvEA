#' @import ggplot2
#'
NULL


#' Runs UMAP on CODEX protein expression. If CleanCODEX has been run,
#' calls UMAP on the cleaned CODEX protein expression. Otherwise,
#' runs UMAP on the original or filtered CODEX protein expression.
#'
#' @param stvea_object STvEA.data class object
#' @param ... parameters to pass into umap
#'
#' @export
#'
GetUmapCODEX <- function(stvea_object, ...) {
  if (!requireNamespace("umap", quietly = TRUE)) {
    stop("Package \"umap\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!is.null(stvea_object@codex_clean)) {
    res <- umap::umap(stvea_object@codex_clean, ...)
    stvea_object@codex_emb <- as.data.frame(res$layout)
    stvea_object@codex_knn <- res$knn$indexes
  } else if (!is.null(stvea_object@codex_protein)) {
    res <- umap::umap(stvea_object@codex_protein, ...)
    stvea_object@codex_emb <- as.data.frame(res$layout)
    stvea_object@codex_knn <- res$knn$indexes
  } else {
    stop("stvea_object does not contain cleaned or raw CODEX protein expression")
  }
  return(stvea_object)
}


#' Runs UMAP on CITE-seq latent space (if provided) or
#' mRNA expression
#'
#' @param stvea_object STvEA.data class object with CITE-seq
#' latent space or mRNA expression
#' @param ... parameters to pass into umap
#'
#' @export
#'
GetUmapCITE <- function(stvea_object, ...) {
  if (!requireNamespace("umap", quietly = TRUE)) {
    stop("Package \"umap\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!is.null(stvea_object@cite_latent)) {
    res <- umap::umap(stvea_object@cite_latent, ...)
    stvea_object@cite_emb <- as.data.frame(res$layout)
  } else if (!is.null(stvea_object@cite_mRNA)) {
    res <- umap::umap(stvea_object@cite_mRNA, ...)
    stvea_object@cite_emb <- as.data.frame(res$layout)
  } else {
    stop("stvea_object does not contain CITE-seq data")
  }
  return(stvea_object)
}


#' Plot expression of a gene or protein in the CITE-seq UMAP or t-SNE embedding
#'
#' @param stvea_object STvEA.data class with CITE-seq expression and embedding
#' @param name gene or protein to plot - can be a character string, index, or vector with up to two names or indices
#' @param type which type of expression data should be plotted? "RNA" or "protein"
#' @param high_color color of high expression on the color ramp
#' @param high_color2 color of high expression for second gene or protein, if provided
#' @param low_color color of low/no expression on the color ramp
#' @param pt_size ggplot geom_point size
#'
#' @param importFrom scales alpha
#' @param importFrom grDevices colorRampPalette
#'
#' @export
#'
PlotExprCITE <- function(stvea_object, name,
                         type="RNA",
                         high_color="red",
                         high_color2 = "green",
                         low_color="light gray",
                         pt_size =0.8) {
  print_type <- type
  if (type == "RNA") {
    plotting_data <- stvea_object@cite_mRNA_norm
  } else if (type == "protein") {
    plotting_data <- stvea_object@cite_protein
    print_type <- "Protein"
  } else {
    stop("type must be either \"RNA\" or \"protein\"", call. =FALSE)
  }

  if (is.null(stvea_object@cite_emb)) {
    stop("No available 2D representation for CITE-seq. Please call GetUmapCITE() or provide other 2D embedding in stvea_object@cite_emb.")
  }
  if (length(name) > 2) {
    stop("name must be at most length 2", call. =FALSE)
  }
  if (!all(name %in% colnames(plotting_data))) {
    stop(sprintf("\"%s\" not in colnames of %s matrix. Change type parameter to select other matrix",
                 paste(name[!name %in% colnames(plotting_data)],collapse="\", \""),type))
  }

  rbPal1 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color,1)), alpha=TRUE)
  color <- rbPal1(100)[as.numeric(cut(plotting_data[,name[1]],breaks = 100))]
  subtitle <- paste("Expression of ", name[1], " (", high_color, ")", sep="")

  if (length(name) == 2) {
    rbPal2 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color2,1)), alpha=TRUE)
    color2 <- rbPal2(100)[as.numeric(cut(plotting_data[,name[2]],breaks = 100))]
    color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m],color2[m]), alpha=TRUE)(3)[2])
    subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, ")", sep="")
  }
  ggplot(stvea_object@cite_emb,
         aes_string(x=colnames(stvea_object@cite_emb)[1],y=colnames(stvea_object@cite_emb)[2],color=factor(1:length(color)))) +
    geom_point(size=pt_size) +
    labs(title = paste(print_type,"expression"), subtitle = subtitle) +
    scale_color_manual(values = alpha(color,1), guide=FALSE) +
    theme_void()
}


#' Plot location of a CITE-seq cell in the CITE-seq UMAP or t-SNE embedding
#'
#' @param stvea_object STvEA.data class with CITE-seq embedding
#' @param index index of CITE-seq cell to plot
#' @param high_color color of cell of interest
#' @param low_color color of all other cells
#' @param pt_size ggplot geom_point size
#'
#' @export
#'
PlotIndexCITE <- function(stvea_object, index,
                          high_color="red",
                          low_color="gray",
                          pt_size =0.8) {
  if (is.null(stvea_object@cite_emb)) {
    stop("No available 2D representation for CITE-seq. Please call GetUmapCITE() or provide other 2D embedding in stvea_object@cite_emb.")
  }
  temp_df <- as.data.frame(rbind(stvea_object@cite_emb[-index,],stvea_object@cite_emb[index,]))
  ggplot(temp_df, aes_string(x=colnames(temp_df)[1],
                             y=colnames(temp_df)[2],
                             color=factor(c(rep(0,nrow(temp_df)-1),1)))) +
    geom_point(size=pt_size) +
    scale_color_manual(values=c(low_color,high_color), guide=FALSE) +
    theme_void()
}


#' Plot clusters of CITE-seq cells in CITE-seq UMAP or t-SNE embedding
#'
#' @param stvea_object STvEA.data class with CITE-seq embedding and clustering
#'
#' @importFrom colorspace rainbow_hcl
#'
#' @export
#'
PlotClusterCITE <- function(stvea_object, pt_size=0.5) {
  if (is.null(stvea_object@cite_emb)) {
    stop("No available 2D representation for CITE-seq. Please call GetUmapCITE() or provide other 2D embedding in stvea_object@cite_emb.")
  }
  if (!length(stvea_object@cite_clusters) == 0) {
    cite_clusters <- stvea_object@cite_clusters
  } else {
    cite_clusters <- rep(-1, nrow(stvea_object@cite_emb))
  }
  if (-1 %in% cite_clusters) {
    colors <- c("gray", rainbow_hcl(length(unique(cite_clusters))-1, c = 80))
  } else {
    colors <- rainbow_hcl(length(unique(cite_clusters)), c = 80)
  }
  ggplot(stvea_object@cite_emb, aes_string(x=colnames(stvea_object@cite_emb)[1],
                                    y=colnames(stvea_object@cite_emb)[2],
                                    color=factor(cite_clusters))) +
    geom_point(size=pt_size) +
    scale_color_manual(values = colors, name="cluster") +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme_void()
}


#' Plot clusters of CODEX cells in CODEX UMAP embedding
#'
#' @param stvea_object STvEA.data class with CODEX embedding and clustering
#'
#' @importFrom colorspace rainbow_hcl
#'
#' @export
#'
PlotClusterCODEXemb <- function(stvea_object, pt_size=0.5) {
  if (is.null(stvea_object@codex_emb)) {
    stop("No available 2D representation for CODEX. Please call GetUmapCODEX() or provide other 2D embedding in stvea_object@codex_emb.")
  }
  if (!length(stvea_object@codex_clusters) == 0) {
    codex_clusters <- stvea_object@codex_clusters
  } else {
    codex_clusters <- rep(-1, nrow(stvea_object@codex_emb))
  }
  if (-1 %in% codex_clusters) {
    colors <- c("gray", rainbow_hcl(length(unique(codex_clusters))-1, c = 80))
  } else {
    colors <- rainbow_hcl(length(unique(codex_clusters)), c = 80)
  }
  ggplot(stvea_object@codex_emb, aes_string(x=colnames(stvea_object@codex_emb)[1],
                                     y=colnames(stvea_object@codex_emb)[2],
                                     color=factor(codex_clusters))) +
    geom_point(size=pt_size) +
    scale_color_manual(values = colors, name="cluster") +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme_void()
}


#' Plot expression of a gene or protein in the CODEX UMAP or t-SNE embedding
#'
#' @param stvea_object STvEA.data class with CODEX expression and embedding
#' @param name gene or protein name to plot
#' @param type which type of expression data should be plotted? "RNA" or "protein"
#' @param high_color color of high expression on the color ramp
#' @param high_color2 color of high expression for the second protein or gene, if provided
#' @param low_color color of low/no expression on the color ramp
#' @param pt_size ggplot geom_point size
#'
#' @export
#'
PlotExprCODEXemb <- function(stvea_object, name,
                              type="protein",
                              high_color="red",
                              high_color2="green",
                              low_color="light gray",
                              pt_size =0.8) {
  print_type <- type
  if (type == "protein") {
    if (!is.null(stvea_object@codex_clean)) {
      plotting_data <- stvea_object@codex_clean
    } else if (!is.null(stvea_object@codex_protein)) {
      plotting_data <- stvea_object@codex_protein
    } else {
      stop("stvea_object must contain CODEX protein data with type=\"protein\"")
    }
    print_type <- "Protein"
  } else if (type == "RNA") {
    plotting_data <- stvea_object@codex_mRNA
  } else {
    stop("type must be either \"RNA\" or \"protein\"", call. =FALSE)
  }

  if (is.null(stvea_object@codex_emb)) {
    stop("No available 2D representation for CODEX. Please call GetUmapCODEX() or provide other 2D embedding in stvea_object@codex_emb.")
  }
  if (length(name) > 2) {
    stop("name must be at most length 2", call. =FALSE)
  }
  if (!all(name %in% colnames(plotting_data))) {
    stop(sprintf("\"%s\" not in colnames of %s matrix. Change type parameter to select other matrix",
                 paste(name[!name %in% colnames(plotting_data)],collapse="\", \""),type))
  }

  rbPal1 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color,1)), alpha=TRUE)
  color <- rbPal1(100)[as.numeric(cut(plotting_data[,name[1]],breaks = 100))]
  subtitle <- paste("Expression of ", name[1], " (", high_color, ")", sep="")

  if (length(name) == 2) {
    rbPal2 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color2,1)), alpha=TRUE)
    color2 <- rbPal2(100)[as.numeric(cut(plotting_data[,name[2]],breaks = 100))]
    color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m],color2[m]), alpha=TRUE)(3)[2])
    subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, ")", sep="")
  }

  ggplot(stvea_object@codex_emb,
         aes_string(x=colnames(stvea_object@codex_emb)[1],y=colnames(stvea_object@codex_emb)[2],
             color=factor(1:length(color)))) +
    geom_point(size=pt_size) +
    labs(title = paste(print_type,"expression"), subtitle = subtitle) +
    scale_color_manual(values = alpha(color,1), guide=FALSE) +
    theme_void()


}


#' Plot expression of a gene or protein in the CODEX spatial coordinates
#'
#' @param stvea_object STvEA.data class with CODEX expression and spatial xy
#' @param name gene or protein name to plot
#' @param type which type of expression data should be plotted? "RNA" or "protein"
#' @param high_color color of high expression on the color ramp
#' @param high_color2 color of high expression for second protein or gene, if provided
#' @param low_color color of low/no expression on the color ramp
#' @param pt_size ggplot geom_point size
#'
#' @export
#'
PlotExprCODEXspatial <- function(stvea_object, name,
                                 type="protein",
                                 high_color="red",
                                 high_color2="green",
                                 low_color="white",
                                 pt_size =0.8) {
  if (is.null(stvea_object@codex_spatial)) {
    stop("stvea_object does not contain CODEX spatial information")
  }
  if (type == "protein") {
    if (!is.null(stvea_object@codex_clean)) {
      plotting_data <- stvea_object@codex_clean
    } else if (!is.null(stvea_object@codex_protein)) {
      plotting_data <- stvea_object@codex_protein
    } else {
      stop("stvea_object must contain CODEX protein data with type=\"protein\"")
    }
  } else if (type == "RNA") {
    plotting_data <- stvea_object@codex_mRNA
  } else {
    stop("type must be either \"RNA\" or \"protein\"", call. =FALSE)
  }

  if (length(name) > 2) {
    stop("name must be at most length 2", call. =FALSE)
  }

  if (!all(name %in% colnames(plotting_data))) {
    stop(sprintf("\"%s\" not in colnames of %s matrix. Change type parameter to select other matrix",
                 paste(name[!name %in% colnames(plotting_data)],collapse="\", \""),type))
  }

  rbPal1 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color,1)), alpha=TRUE)
  color <- rbPal1(100)[as.numeric(cut(plotting_data[,name[1]],breaks = 100))]
  subtitle <- paste("Expression of ", name[1], " (", high_color, ")", sep="")

  if (length(name) == 2) {
    rbPal2 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color2,1)), alpha=TRUE)
    color2 <- rbPal2(100)[as.numeric(cut(plotting_data[,name[2]],breaks = 100))]
    color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m],color2[m]), alpha=TRUE)(3)[2])
    subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, ")", sep="")
  }

  x_tmp <- stvea_object@codex_spatial[,"x"]
  x_tmp <- x_tmp - min(x_tmp)
  y_tmp <- stvea_object@codex_spatial[,"y"]
  y_tmp <- y_tmp - min(y_tmp)
  spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp))
  ggplot(spatial_tmp,
         aes(x=x,y=y,color=factor(1:length(color)))) +
    geom_point(size=0.5, alpha=0.5) +
    scale_color_manual(values=alpha(color,1)) +
    guides(color=FALSE) +
    ylim(max(y_tmp), 0) +
    labs(title = paste("Spatial",type,"expression"), subtitle = subtitle) +
    theme_void() + coord_fixed()
}


#' Plot location of a CODEX cell in the CODEX UMAP or t-SNE embedding
#'
#' @param stvea_object STvEA.data class with CODEX embedding
#' @param index index of CODEX cell to plot
#' @param high_color color of cell of interest
#' @param low_color color of all other cells
#' @param pt_size ggplot geom_point size
#'
#' @export
#'
PlotIndexCODEXemb <- function(stvea_object, index,
                               high_color="red",
                               low_color="gray",
                               pt_size =0.8) {
  if (is.null(stvea_object@codex_emb)) {
    stop("No available 2D representation for CODEX. Please call GetUmapCODEX() or provide other 2D embedding in stvea_object@codex_emb.")
  }
  temp_df <- as.data.frame(rbind(stvea_object@codex_emb[-index,],stvea_object@codex_emb[index,]))
  ggplot(temp_df,
         aes_string(x=colnames(temp_df)[1],y=colnames(temp_df)[2],color=factor(c(rep(0,nrow(temp_df)-1),1)))) +
    geom_point(size=pt_size) +
    scale_color_manual(values=c(low_color,high_color), guide=FALSE) +
    theme_void()
}


#' Plot location of a CODEX cell in the CODEX spatial coordinates
#'
#' @param stvea_object STvEA.data class with CODEX spatial xy
#' @param index index of CODEX cell to plot
#' @param high_color color of cell of interest
#' @param low_color color of all other cells
#' @param pt_size ggplot geom_point size
#'
#' @export
#'
PlotIndexCODEXspatial <- function(stvea_object, index,
                                  high_color="red",
                                  low_color="gray",
                                  pt_size =0.8) {
  if (is.null(stvea_object@codex_spatial)) {
    stop("stvea_object does not contain CODEX spatial information")
  }
  temp_df <- as.data.frame(rbind(stvea_object@codex_spatial[-index,],stvea_object@codex_spatial[index,]))
  ggplot(temp_df,
         aes(x=x,y=y,color=c(rep("other",nrow(temp_df)-1),"this"))) +
    geom_point(size=pt_size) +
    scale_color_manual(values=c(low_color,high_color), guide=FALSE) +
    guides(color = FALSE) +
    ylim(max(temp_df$y), 0) +
    theme_void() + coord_fixed()
}
