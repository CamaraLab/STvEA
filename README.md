# STvEA

```STvEA``` is an analysis pipeline for cleaning, clustering, and plotting CODEX and CITE-seq protein data, mapping a CODEX dataset to a matching CITE-seq dataset, and assessing colocalization of features using the Adjacency Score. More information can be found in:

K. W. Govek*, E. C. Troisi*, Z. Miao, R. G. Aubin, S. Woodhouse, and P. G. Cámara. _Single-Cell Transcriptomic Analysis of mIHC Images via Antigen Mapping_. **Science Advances** 7 (2021) 10. [DOI: 10.1126/sciadv.abc5464](https://doi.org/10.1126/sciadv.abc5464). *authors contributed equally.

## Installation
Install flowCore
```
# Bioconductor 3.6 (R 3.4)
source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")

------ or ------

# Bioconductor 3.8+ (R 3.5+)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("flowCore")
```
Install STvEA
```
devtools::install_github("CamaraLab/STvEA")
```
## Docker image
We provide Docker images to run STvEA in R Studio, based off the rocker/rstudio images:
https://hub.docker.com/r/camaralab/stvea

## Tutorials

[Analyzing CODEX data](https://github.com/CamaraLab/STvEA/tree/master/examples/codex_tutorial.md)

Example pipeline for only CODEX data. Performs filtering of protein expression data, visualization of clusters, and analysis of the spatial co-localization of both clusters and proteins.

[Analyzing CyTOF data and reading from FCS files](https://github.com/CamaraLab/STvEA/tree/master/examples/cytof_tutorial.md)

Example pipeline for reading CyTOF or CODEX data from FCS files and analyzing mIHC or cytometry data without a spatial component.

[Mapping CODEX data to CITE-seq data](https://github.com/CamaraLab/STvEA/tree/master/examples/mapping_tutorial.md)

Performs separate filtering and clustering for both CODEX and CITE-seq protein expression data. Maps the CODEX protein space to the CITE-seq protein space in order to transfer features such as gene expression and clusters from the CITE-seq cells onto the CODEX cells. Analyzes co-localization of genes, proteins, and clusters using the Adjacency Score.

[Mapping CODEX data to CITE-seq data analyzed using Seurat](https://github.com/CamaraLab/STvEA/tree/master/examples/seurat_tutorial.md)

Retrieves information about the CITE-seq mRNA and protein expression data from a Seurat object. Performs separate filtering and clustering for both CODEX and CITE-seq protein expression data. Maps the CODEX protein space to the CITE-seq protein space in order to transfer features such as gene expression and clusters from the CITE-seq cells onto the CODEX cells. Analyzes co-localization of genes, proteins, and clusters using the Adjacency Score.
