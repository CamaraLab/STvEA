#' @import ggplot2
#'
NULL

#' Plot expression of a gene or protein in the CITE-seq UMAP embedding
#'
#' @param stvea_object STvEA.data class with CITE-seq expression and embedding
#' @param name gene or protein name to plot
#' @param type which type of expression data should be plotted? "RNA" or "protein"
#'
#' @export
#'
PlotExprCITE <- function(stvea_object, name, type="RNA") {
  if (type == "RNA") {
    ggplot(stvea_object@cite_emb,
           aes(x=V1,y=V2,color=stvea_object@cite_mRNA[,name])) +
      geom_point(size=0.8) +
      scale_color_gradient(low="gray", high="red", guide=FALSE) +
      theme_void()
  } else if (type == "protein") {
    ggplot(stvea_object@cite_emb,
           aes(x=V1,y=V2,color=stvea_object@cite_protein[,name])) +
      geom_point(size=0.8) +
      scale_color_gradient(low="gray", high="blue", guide=FALSE) +
      theme_void()
  } else {
    stop("Invalid type setting")
  }
}


#' Plot location of a CITE-seq cell in the CITE-seq UMAP embedding
#'
#' @param stvea_object STvEA.data class with CITE-seq embedding
#' @param index index of CITE-seq cell to plot
#'
#' @export
#'
PlotIndexCITE <- function(stvea_object, index) {
  temp_df <- as.data.frame(rbind(stvea_object@cite_emb[-index,],stvea_object@cite_emb[index,]))
  ggplot(temp_df,
         aes(x=V1,y=V2,color=c(rep("other",nrow(temp_df)-1),"this"))) +
    geom_point(size=0.8) +
    scale_color_manual(values=c("gray","red"), guide=FALSE) +
    theme_void()
}


#' Plot expression of a gene or protein in the CODEX UMAP embedding
#'
#' @param stvea_object STvEA.data class with CODEX expression and embedding
#' @param name gene or protein name to plot
#' @param type which type of expression data should be plotted? "RNA" or "protein"
#'
#' @export
#'
PlotExprCODEXumap <- function(stvea_object, name, type="protein") {
  if (type == "RNA") {
    ggplot(stvea_object@codex_emb,
           aes(x=V1,y=V2,color=stvea_object@codex_mRNA[,name])) +
      geom_point(size=0.8) +
      scale_color_gradient(low="gray", high="red", guide=FALSE) +
      theme_void()
  } else if (type == "protein") {
    ggplot(stvea_object@codex_emb,
           aes(x=V1,y=V2,color=stvea_object@codex_protein[,name])) +
      geom_point(size=0.8) +
      scale_color_gradient(low="gray", high="purple", guide=FALSE) +
      theme_void()
  } else {
    stop("Invalid type setting")
  }
}


#' Plot expression of a gene or protein in the CODEX spatial coordinates
#'
#' @param stvea_object STvEA.data class with CODEX expression and spatial xy
#' @param name gene or protein name to plot
#' @param type which type of expression data should be plotted? "RNA" or "protein"
#'
#' @export
#'
PlotExprCODEXspatial <- function(stvea_object, name, type="protein") {
  if (type == "protein") {
    ggplot(as.data.frame(stvea_object@codex_spatial),
           aes(x=x,y=y,color=stvea_object@codex_protein[,name])) +
      geom_point(size=0.5) +
      scale_color_gradient(low="white", high="blue") +
      guides(color = FALSE) +
      ylim(max(y), 0) +
      theme_void()
  } else if (type == "RNA") {
    ggplot(as.data.frame(stvea_object@codex_spatial),
           aes(x=x,y=y,color=stvea_object@codex_mRNA[,name])) +
      geom_point(size=0.5) +
      scale_color_gradient(low="white", high="blue") +
      guides(color = FALSE) +
      ylim(max(y), 0) +
      theme_void()
  } else {
    stop("Invalid type setting")
  }
}


#' Plot location of a CODEX cell in the CODEX UMAP embedding
#'
#' @param stvea_object STvEA.data class with CODEX embedding
#' @param index index of CODEX cell to plot
#'
#' @export
#'
PlotIndexCODEXumap <- function(stvea_object, index) {
  temp_df <- as.data.frame(rbind(stvea_object@codex_emb[-index,],stvea_object@codex_emb[index,]))
  ggplot(temp_df,
         aes(x=V1,y=V2,color=c(rep("other",nrow(temp_df)-1),"this"))) +
    geom_point(size=0.8) +
    scale_color_manual(values=c("gray","red"), guide=FALSE) +
    theme_void()
}


#' Plot location of a CODEX cell in the CODEX spatial coordinates
#'
#' @param stvea_object STvEA.data class with CODEX spatial xy
#' @param index index of CODEX cell to plot
#'
#' @export
#'
PlotIndexCODEXspatial <- function(stvea_object, index) {
  temp_df <- as.data.frame(rbind(stvea_object@codex_spatial[-index,],stvea_object@codex_spatial[index,]))
  ggplot(temp_df,
         aes(x=x,y=y,color=c(rep("other",nrow(temp_df)-1),"this"))) +
    geom_point(size=0.5) +
    scale_color_gradient(low="gray", high="red") +
    guides(color = FALSE) +
    ylim(max(y), 0) +
    theme_void()

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
