
#'
#'
# mapping_data_holder <- setClass(
#   Class = 'mapping_data_holder',
#   slots = c(
#     cite_latent = 'ANY', # cite x low
#     cite_protein = 'ANY', # cite x 30
#     cite_protein_clean = 'ANY', # cite x 30
#
#     codex_protein = 'ANY', # codex x 30
#     codex_protein_filter = 'ANY', # codex x 30
#     codex_protein_clean = 'ANY', # codex x 30
#
#     corrected_codex = 'ANY', # codex x 30 (can only do this if CITE always ref)
#     codex_nn = 'dgCMatrix', # cite x codex
#     cite_nn = 'dgCMatrix', # codex x cite
#   )
# )


# other matrices: cite_gene, codex_spatial, codex_gene,
# cite_cluster_labels, codex_cluster_labels, cite_umap_embedding

# maybe have easy method for transfering data between holders?
