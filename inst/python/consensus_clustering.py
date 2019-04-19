import hdbscan
import numpy as np
from sklearn.metrics import silhouette_score

def run_hdbscan(latent, umap_latent, min_cluster_size_list, min_sample_list, cache_dir=None):
    latent = np.array(latent)
    num_cells = latent.shape[0]
    umap_latent = np.array(umap_latent)
    #min_cluster_size_list = list(range(10,40,6))
    #min_sample_list = list(range(30,80,5))

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

