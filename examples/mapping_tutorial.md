Using STvEA to map CODEX and CITE-seq data
================

``` r
library(STvEA)
```

    ## Warning: replacing previous import 'Matrix::expm' by 'expm::expm' when
    ## loading 'AdjacencyScore'

Read in CODEX data
------------------

Dataframe converted from FCS files at <http://welikesharingdata.blob.core.windows.net/forshare/index.html>

``` r
data("codex_balbc1")
```

### Take corner section of CODEX data

The clustering and Adjacency Score functions are fairly slow on large numbers of cells.

``` r
codex_subset <- codex_spatial$x < 3000 & codex_spatial$y < 3000
codex_protein <- codex_protein[codex_subset,]
codex_blanks <- codex_blanks[codex_subset,]
codex_size <- codex_size[codex_subset]
codex_spatial <- codex_spatial[codex_subset,]
```

Read in CITE-seq data
---------------------

``` r
# link to CITE-seq data on Dropbox
```

Create object to hold data
--------------------------

The STvEA.data class conveniently handles the required data frames and matrices between function calls.

``` r
stvea_object <- SetDataCODEX(codex_protein = codex_protein,
                             codex_blanks = codex_blanks,
                             codex_size = codex_size,
                             codex_spatial = codex_spatial)
stvea_object <- SetDataCITE(cite_mRNA = cite_mRNA,
                            cite_protein = cite_protein,
                            cite_latent = cite_latent,
                            stvea_object = stvea_object)
```

Filter and clean CODEX protein expression
-----------------------------------------

Remove cells that are too small or large, or have too low or too high expression in the blank channels. If lower and upper limits aren't specified, quantiles are taken as the limits. Then normalize data by the total counts per cell

``` r
stvea_object <- FilterCODEX(stvea_object, size_lim = c(1000, 25000),
                            blank_lower = c(-1200, -1200, -1200, -1200),
                            blank_upper = c(6000, 2500, 5000, 2500))
```

Fit a Gaussian mixture model to the expression levels of each protein. New data will be the cumulative probability according to the Gaussian with the higher mean

``` r
stvea_object <- CleanCODEX(stvea_object)
```

Clean and normalize CITE-seq protein expression
-----------------------------------------------

Fit a Negative Binomial mixture model to the expression levels of each protein in the CITE-seq dataset. Cleaned data will be the cumulative probability according to the Negative Binomial with the higher median.

``` r
# This will take around 10 minutes
stvea_object <- CleanCITE(stvea_object, num_cores=8)
```

Cluster CITE-seq cells based on mRNA expression
-----------------------------------------------

Compute the 2 dimensional UMAP embedding of the CITE-seq mRNA latent space for later visualization.

``` r
stvea_object <- GetUmapCITE(stvea_object, n_neighbors = 50, min_dist=0.1, negative_sample_rate=50)
```

Run UMAP on the CITE-seq latent space to the same number of dimensions (running UMAP before density based clustering is suggested in the UMAP documentation <https://umap-learn.readthedocs.io/en/latest/clustering.html>). Scan over the two parameters of the Python HDBSCAN implementation, min\_cluster\_size and min\_samples (more information on selecting these paramters <https://hdbscan.readthedocs.io/en/latest/parameter_selection.html>).

``` r
stvea_object <- ParameterScan(stvea_object, min_cluster_size_range = seq(5,20,4), min_sample_range = seq(10,40,3), n_neighbors=50, min_dist=0.1, negative_sample_rate=50)
```

    ## Running UMAP on the CITE-seq latent space
    ## Running HDBSCAN on the UMAP space

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-11-1.png)

Create a dissimilarity matrix between cells as the number of HDBSCAN clusterings (passing a certain silhouette cutoff) in which they are assigned to the same cluster. Perform agglomerative clustering on this dissimilarity matrix using Python's scipy.cluster.hierarchy, cutting the tree using the inconsistent value. Remove all clusters that have fewer than 10 cells.

``` r
stvea_object <- ConsensusCluster(stvea_object, silhouette_cutoff = 0.1, inconsistent_value = 0.3, min_cluster_size = 10)
```

Cluster CODEX cells based on protein expression
-----------------------------------------------

Compute the 2 dimensional UMAP embedding of the cleaned CODEX protein expression for later visualization. UMAP also returns the KNN indices with k = n\_neighbors.

``` r
# This will take around 5 minutes for ~10000 cells
stvea_object <- GetUmapCODEX(stvea_object, metric = 'pearson', n_neighbors=30,
                             min_dist=0.1, negative_sample_rate = 50)
```

Use the KNN indices from UMAP to perform Louvain clustering

``` r
stvea_object <- ClusterCODEX(stvea_object, k=30)
```

Map CITE-seq and CODEX protein
------------------------------

Perform a modified version of Seurat anchor correction from <https://www.biorxiv.org/content/10.1101/460147v1> mapping the CODEX protein space into the CITE-seq protein space.

``` r
stvea_object <- MapCODEXtoCITE(stvea_object, num_chunks=8, seed=30, num_cores=4)
```

Create neighbor matrices
------------------------

Create transfer matrices of k CITE-seq nearest neighbors for each CODEX cell and k CODEX nearest neighbors for each CITE-seq cell after the anchor correction above

``` r
stvea_object <- GetTransferMatrix(stvea_object)
```

Visualize clustering, mRNA, and protein expression
--------------------------------------------------

Color each cell in the CITE-seq UMAP embedding with its cluster assignment. Cells in gray were not assigned to any cluster.

``` r
PlotClusterCITE(stvea_object)
```

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-17-1.png)

Color each cell in the CITE-seq UMAP embedding with its expression level of one of two genes. If two gene names are provided, color will be interpolated between red and green color values.

``` r
PlotExprCITE(stvea_object, c("Cd4", "Ighd"), type="RNA")
```

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-18-1.png)

Color each cell in the CITE-seq UMAP embedding with its expression level of one of two proteins.

``` r
PlotExprCITE(stvea_object, c("CD4","IgD"), type="protein")
```

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-19-1.png)

Color each cell in the CODEX UMAP embedding with its cluster assignment. Cells in gray were not assigned to any cluster.

``` r
PlotClusterCODEXemb(stvea_object)
```

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-20-1.png)

Color each cell in the CODEX UMAP embedding with its expression level of one or two proteins.

``` r
PlotExprCODEXemb(stvea_object, c("CD4","IgD"))
```

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-21-1.png)

Color the CODEX spatial slide with the expression level of one or two proteins.

``` r
PlotExprCODEXspatial(stvea_object, c("CD4","IgD"))
```

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-22-1.png)

Color the CODEX spatial slide with the expression level of one or two genes that were mapped from the CITE-seq expression levels.

``` r
PlotExprCODEXspatial(stvea_object, c("Cd4", "Ighd"), type="RNA")
```

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-23-1.png)

Assess colocalization of clusters
---------------------------------

Run the Adjacency Score (<https://github.com/CamaraLab/AdjacencyScore>) to evaluate how often two features take high values in adjacent nodes in a KNN graph of the CODEX spatial dimensions.

Assess which pairs of CODEX clusters often appear in adjacent cells

``` r
codex_cluster_adj <- AdjScoreClustersCODEX(stvea_object, k=3)
```

    ## Creating permutation matrices - 0.007 seconds
    ## Computing adjacency score for each feature pair - 0.389 seconds

``` r
AdjScoreHeatmap(codex_cluster_adj)
```

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-24-1.png)

Assess which pairs of clusters from the CITE-seq mRNA analysis are often mapped to neighboring CODEX cells

``` r
cite_cluster_adj <- AdjScoreClustersCITE(stvea_object, k=3, num_cores=8)
```

    ## Creating permutation matrices - 7.001 seconds
    ## Computing adjacency score for each feature pair - 31.34 seconds

``` r
AdjScoreHeatmap(cite_cluster_adj)
```

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-25-1.png)

Assess which pairs of proteins are often highly expressed in adjacent cells

``` r
protein_adj <- AdjScoreProteins(stvea_object, k=3, num_cores=8)
```

    ## Creating permutation matrices - 8.971 seconds
    ## Computing adjacency score for each feature pair - 37.335 seconds

``` r
AdjScoreHeatmap(protein_adj)
```

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-26-1.png)

Assess which pairs of genes from the given list are often highly expressed in adjacency cells. There are too many genes to compute the Adjacency Score for all combinations, so let's compute it for all pairs of genes corresponding to one of the proteins.

``` r
gene_list <- c("Ptprc", "Ly6c1", "Trdc", "Ly6g", "Cd19", "Vcam1", "Cd3g",
               "Fcgr3", "Cd8a", "Thy1", "Adgre1", "Itgax", "Car1", "Itgam",
               "Ighd", "Cd27", "Cd5", "Cd79b", "Tfrc", "Pecam1", "Cd4",
               "Ighm", "Cr2", "Cd44", "Ncr1", "H2-Ab1")
gene_pairs <- t(combn(gene_list,2))
for (gene in gene_list) {
  gene_pairs <- rbind(gene_pairs, c(gene,gene))
}
gene_adj <- AdjScoreGenes(stvea_object, gene_pairs,  k=3, num_cores=8)
```

    ## Creating permutation matrices - 9.451 seconds
    ## Computing adjacency score for each feature pair - 37.57 seconds

``` r
AdjScoreHeatmap(gene_adj)
```

![](mapping_tutorial_files/figure-markdown_github/unnamed-chunk-27-1.png)