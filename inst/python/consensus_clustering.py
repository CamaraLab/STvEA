import hdbscan
import numpy as np
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

def run_hdbscan(latent, umap_latent, min_cluster_size_list, min_sample_list, cache_dir=None):
    latent = np.array(latent)
    num_cells = latent.shape[0]
    umap_latent = np.array(umap_latent)

    results = []
    for min_samples in min_sample_list:
        for min_cluster_size in min_cluster_size_list:
            if (cache_dir is None):
                clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples, metric="correlation")
            else:
                clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples, metric="correlation", memory=cache_dir)
            hdbscan_labels = clusterer.fit_predict(umap_latent)
            score = silhouette_score(latent, hdbscan_labels)
            results.append(hdbscan_labels.tolist())
    return results

def consensus_cluster(distance_matrix, inconsistent_value = 0.3, min_cluster_size=10):
    new_distance = squareform(distance_matrix)
    hierarchical_tree = linkage(new_distance, "average")
    hier_consensus_labels = fcluster(hierarchical_tree, t=inconsistent_value)
    hier_unique_labels = set(hier_consensus_labels)
    for label in hier_unique_labels:
        indices = [i for i, x in enumerate(hier_consensus_labels) if x == label]
        if len(indices) < min_cluster_size:
            for index in indices:
                hier_consensus_labels[index] = -1
    print(type(hier_consensus_labels))
    return hier_consensus_labels.tolist()
