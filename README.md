# STvEA

```STvEA``` is an analysis pipeline for cleaning, clustering, and plotting CODEX and CITE-seq protein data, mapping a CODEX dataset to a matching CITE-seq dataset, and assessing colocalization of features using the Adjacency Score. More information can be found in:

- K. W. Govek*, E. C. Troisi*, Z. Miao, R. G. Aubin, S. Woodhouse, and P. G. Cámara. _Single-Cell Transcriptomic Analysis of mIHC Images via Antigen Mapping._ Science Advances (2021), **In press**. [DOI: 10.1101/672501 (bioRxiv)](https://www.biorxiv.org/content/10.1101/672501v2). *authors contributed equally.

## Installation

```
devtools::install_github("CamaraLab/STvEA")
```

## Tutorials

[Analyzing CODEX data](https://github.com/CamaraLab/STvEA/tree/master/examples/codex_tutorial.md)

Example pipeline for only CODEX data. Performs filtering of protein expression data, visualization of clusters, and analysis of the spatial co-localization of both clusters and proteins.

[Mapping CODEX data to CITE-seq data](https://github.com/CamaraLab/STvEA/tree/master/examples/mapping_tutorial.md)

Performs separate filtering and clustering for both CODEX and CITE-seq protein expression data. Maps the CODEX protein space to the CITE-seq protein space in order to transfer features such as gene expression and clusters from the CITE-seq cells onto the CODEX cells. Analyzes co-localization of genes, proteins, and clusters using the Adjacency Score.

[Mapping CODEX data to CITE-seq data analyzed using Seurat](https://github.com/CamaraLab/STvEA/tree/master/examples/seurat_tutorial.md)

Retrieves information about the CITE-seq mRNA and protein expression data from a Seurat object. Performs separate filtering and clustering for both CODEX and CITE-seq protein expression data. Maps the CODEX protein space to the CITE-seq protein space in order to transfer features such as gene expression and clusters from the CITE-seq cells onto the CODEX cells. Analyzes co-localization of genes, proteins, and clusters using the Adjacency Score.
