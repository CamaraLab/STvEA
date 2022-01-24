library(purrr)
library(Matrix)

cite_umap_emb = read.table("umap_embedding_v2.csv", sep=",", header=FALSE, row.names=1, stringsAsFactors = FALSE)
cite_protein <- read.csv("protein_expr_corrected_clean.csv", row.names=1, stringsAsFactors = FALSE) %>% as.matrix
colnames(cite_protein) <- colnames(cite_protein) %>% map_chr(function(x) sub("ADT_", "", x))

spatial1 = read.table("CODEX_BALBc-1_filtered_spatial.txt")
spatial2 = read.table("CODEX_BALBc-2_filtered_spatial.txt")
spatial3 = read.table("CODEX_BALBc-3_filtered_spatial.txt")

codex_exp_clean_BALBc1 = readRDS("codex_exp_clean_BALBc1.rds") %>% as.matrix
codex_exp_clean_BALBc2 = readRDS("codex_exp_clean_BALBc2.rds") %>% as.matrix
codex_exp_clean_BALBc3 = readRDS("codex_exp_clean_BALBc3.rds") %>% as.matrix

cite_nn1 = readRDS("cite_nn_exp_BALBc-1.rds")
cite_nn2 = readRDS("cite_nn_exp_BALBc-2.rds")
cite_nn3 = readRDS("cite_nn_exp_BALBc-3.rds")

cite_gene = readRDS("cite_gene.rds") # normalised, but not log transformed
# genes = read.table("genes.tsv", row.names=1, header=F, sep="\t", stringsAsFactors=F)
## cite_gene <- read.csv("gene_matrix_all_new.csv",row.names=1, stringsAsFactors=F) %>% as.matrix
# cite_gene <- read.csv("../../../Data/CITE-seq/balbc_all/gene_matrix_all_uncut.csv",row.names=1, stringsAsFactors=F) %>% as.matrix
# cite_gene <- t(t(cite_gene)/colSums(cite_gene))
# row_names <- genes[row.names(cite_gene),1]
# cite_gene <- cite_gene[!duplicated(row_names),]
# row_names <- row_names[!duplicated(row_names)]
# row.names(cite_gene) <- row_names
# cite_gene = Matrix(cite_gene, sparse=T)
# saveRDS(cite_gene, "cite_gene.rds")

clusters = read.csv("cluster_labels.csv", header=F, row.names=1)

guide.image1 = png::readPNG("www/codex-mouse1.png")
guide.image2 = png::readPNG("www/codex-mouse2.png")
guide.image3 = png::readPNG("www/codex-mouse3.png")

# for pilot data
# genes = read.table("genes.tsv", sep="\t", header=FALSE, stringsAsFactors = FALSE, row.names=1)
# tcounts <- read.table("../../../Data/CITE-seq/black6/filtered_gene_matrix_combined_uncut.csv", sep=",", header=TRUE, stringsAsFactors = FALSE, row.names=1) %>% as.matrix
# row_names <- genes[row.names(tcounts),1]
# tcounts <- tcounts[!duplicated(row_names),]
# row_names <- row_names[!duplicated(row_names)]
# row.names(tcounts) <- row_names
# saveRDS(tcounts, "filtered_gene_matrix_combined_uncut.rds")
tcounts <- readRDS("filtered_gene_matrix_combined_uncut.rds")
tumap = read.csv("umap_embedding_combined_nocycle.csv", header = FALSE, row.names = 1)
