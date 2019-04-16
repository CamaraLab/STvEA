import ctypes
import hdbscan
import itertools
from multiprocessing import RawArray, Pool
import numpy as np
from scipy.sparse import csr_matrix
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
            rows = []
            cols = []
            data = []
            for cell1 in range(num_cells):
                for cell2 in range(num_cells):
                    if cell1 == cell2:
                        rows.append(cell1)
                        cols.append(cell2)
                        data.append(-1)
                    elif hdbscan_labels[cell1] != -1 and hdbscan_labels[cell1] == hdbscan_labels[cell2]:
                        rows.append(cell1)
                        cols.append(cell2)
                        data.append(-1)
            return_dict = {"rows" : rows, "cols": cols, "data":data, "score":score, "num_cells":num_cells}
            results.append(return_dict)
        
    return results
        
