#' @import ggplot2
#'
NULL

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
  if (type == "RNA") {
    plotting_data <- stvea_object@cite_mRNA_norm
  } else if (type == "protein") {
    plotting_data <- stvea_object@cite_protein
  } else {
    stop("type must be either \"RNA\" or \"protein\"", call. =FALSE)
  }

  if (length(name) > 2) {
    stop("name must be at most length 2", call. =FALSE)
  }

  rbPal1 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color,1)), alpha=TRUE)
  color <- rbPal1(100)[as.numeric(cut(plotting_data[,name[1]],breaks = 100))]

  if (length(name) == 2) {
    rbPal2 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color2,1)), alpha=TRUE)
    color2 <- rbPal2(100)[as.numeric(cut(plotting_data[,name[2]],breaks = 100))]
    color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m],color2[m]), alpha=TRUE)(3)[2])
  }
  ggplot(stvea_object@cite_emb,
         aes_string(x=colnames(stvea_object@cite_emb)[1],y=colnames(stvea_object@cite_emb)[2],color=factor(1:length(color)))) +
    geom_point(size=pt_size) +
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
PlotClusterCITE <- function(stvea_object) {
  if (-1 %in% stvea_object@cite_clusters) {
    colors <- c("gray", rainbow_hcl(length(unique(stvea_object@cite_clusters))-1, c = 80))
  } else {
    colors <- rainbow_hcl(length(unique(stvea_object@cite_clusters)), c = 80)
  }
  ggplot(stvea_object@cite_emb, aes(x=V1,y=V2,color=factor(stvea_object@cite_clusters))) +
    geom_point(size=0.5) +
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
  if (type == "protein") {
    plotting_data <- stvea_object@codex_protein
  } else if (type == "RNA") {
    plotting_data <- stvea_object@codex_mRNA
  } else {
    stop("type must be either \"RNA\" or \"protein\"", call. =FALSE)
  }

  if (length(name) > 2) {
    stop("name must be at most length 2", call. =FALSE)
  }

  rbPal1 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color,1)), alpha=TRUE)
  color <- rbPal1(100)[as.numeric(cut(plotting_data[,name[1]],breaks = 100))]

  if (length(name) == 2) {
    rbPal2 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color2,1)), alpha=TRUE)
    color2 <- rbPal2(100)[as.numeric(cut(plotting_data[,name[2]],breaks = 100))]
    color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m],color2[m]), alpha=TRUE)(3)[2])
  }
  ggplot(stvea_object@codex_emb,
         aes_string(x=colnames(stvea_object@codex_emb)[1],y=colnames(stvea_object@codex_emb)[2],
             color=factor(1:length(color)))) +
    geom_point(size=pt_size) +
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
  if (type == "protein") {
    plotting_data <- stvea_object@codex_protein
  } else if (type == "RNA") {
    plotting_data <- stvea_object@codex_mRNA
  } else {
    stop("type must be either \"RNA\" or \"protein\"", call. =FALSE)
  }

  if (length(name) > 2) {
    stop("name must be at most length 2", call. =FALSE)
  }

  rbPal1 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color,1)), alpha=TRUE)
  color <- rbPal1(100)[as.numeric(cut(plotting_data[,name[1]],breaks = 100))]

  if (length(name) == 2) {
    rbPal2 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color2,1)), alpha=TRUE)
    color2 <- rbPal2(100)[as.numeric(cut(plotting_data[,name[2]],breaks = 100))]
    color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m],color2[m]), alpha=TRUE)(3)[2])
  }

  ggplot(as.data.frame(spatial),
         aes(x=x,y=y,color=factor(1:length(color)))) +
    geom_point(size=0.5, alpha=0.5) +
    scale_color_manual(values=alpha(color,1)) +
    guides(color = FALSE) +
    ylim(max(y), 0) +
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
  temp_df <- as.data.frame(rbind(stvea_object@codex_spatial[-index,],stvea_object@codex_spatial[index,]))
  ggplot(temp_df,
         aes(x=x,y=y,color=c(rep("other",nrow(temp_df)-1),"this"))) +
    geom_point(size=pt_size) +
    scale_color_manual(values=c(low_color,high_color), guide=FALSE) +
    guides(color = FALSE) +
    ylim(max(y), 0) +
    theme_void() + coord_fixed()
}


#' Plot the CODEX neighbors of a set of CITE-seq cells
#' in the spatial coordinates
#'
PlotSpatialNeighbors <- function(cite_cells) {
}

#' Plot CITE-seq neighbors of a set of CODEX cells
#' in the CITE-seq mRNA UMAP
#'
PlotUmapNeighbors <- function(codex_cells) {
}

#' Plot CITE-seq gene expression on
#' CODEX spatial coordinates
#'
PlotSpatialRNA <- function(gene_name) {
}

#' Plot heatmap of cluster feature scores
#'
plot_heatmap <- function(feature_score) {
}
