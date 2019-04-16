import ctypes
import hdbscan
import itertools
from multiprocessing import RawArray, Pool
import numpy as np
from scipy.sparse import csr_matrix
from sklearn.metrics import silhouette_score


def single_clustering_test(arguments):
    print("single test")
    min_cluster_size, min_samples = arguments
    latent = np.frombuffer(var_dict['latent']).reshape(var_dict['latent_shape'])
    umap_latent = np.frombuffer(var_dict['umap_latent']).reshape(var_dict['umap_shape'])
    num_cells = var_dict['latent_shape'][0]
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples, metric="correlation")
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
    distance_matrix = csr_matrix((data,(rows,cols)), shape=(num_cells,num_cells))
    return distance_matrix, score
    
def test_fn(arguments):
    print("in test fn")
    return arguments
    
def init_worker(latent, umap_latent):
    global var_dict

    latent_shared = RawArray('d', latent.shape[0]*latent.shape[1])
    umap_shared = RawArray('d', umap_latent.shape[0]*umap_latent.shape[1])
    latent_np = np.frombuffer(latent_shared, dtype=np.float64).reshape(latent.shape)
    umap_np = np.frombuffer(umap_shared, dtype=np.float64).reshape(umap_latent.shape)
    np.copyto(latent_np, latent)
    np.copyto(umap_np, umap_latent)
    
    var_dict = {}
    var_dict['latent'] = latent
    var_dict['latent_shape'] = latent.shape
    var_dict['umap_latent'] = umap_latent
    var_dict['umap_shape'] = umap_latent.shape
    print("done with init")

def run_hdbscan(latent, umap_latent, min_cluster_size_list, min_sample_list):
    latent = np.array(latent)
    umap_latent = np.array(umap_latent)
    #min_cluster_size_list = list(range(10,40,6))
    #min_sample_list = list(range(30,80,5))
    arguments = list(itertools.product(min_cluster_size_list, min_sample_list))
    
    with Pool(processes = 1, initializer=init_worker, initargs=(latent, umap_latent)) as pool:
        print("start parallel")
        parallel_results = pool.map(test_fn, arguments)
        print(type(parallel_results))
        pool.close()
        pool.join()
        
    return parallel_results

        
def silhouette_cut(parallel_results, silhouette_score_cutoff):
    num_cells = latent.shape[0]
    distance_matrix = csr_matrix((num_cells,num_cells))
    all_scores = []
    total_runs = 0.0
    for result in parallel_results:
        if result[1] > silhouette_score_cutoff:
            distance_matrix = distance_matrix + result[0]
            all_scores.append(result[1])
            total_runs += 1
    distance_matrix = distance_matrix.toarray()
    distance_matrix = distance_matrix + total_runs
    if total_runs > 0:
        distance_matrix /= total_runs
    else:
        print("0 runs selected")
    return distance_matrix, all_scores
        
        
def consensus_cluster(distance_matrix, inconsistent = 0.5, min_cluster_size = 10):
    new_distance = squareform(distance_matrix)
    hierarchical_tree = linkage(new_distance, "average")
    hier_consensus_labels = fcluster(hierarchical_tree, t=inconsistent)
    hier_unique_labels = set(hier_consensus_labels)
    for label in hier_unique_labels:
        indices = [i for i, x in enumerate(hier_consensus_labels) if x == label]
        if len(indices) < min_cluster_size:
            for index in indices:
                hier_consensus_labels[index] = -1
    return hier_consensus_labels
            
