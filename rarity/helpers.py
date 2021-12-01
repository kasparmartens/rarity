import numpy as np

def summarise_binary_profiles(binary_profiles):
    unique_profiles, cluster_allocations, counts = np.unique(binary_profiles, axis=0, return_inverse=True, return_counts=True)
    reorder_inds = counts.argsort()[::-1]
    cluster_allocations2 = np.zeros_like(cluster_allocations)
    for i, e in enumerate(reorder_inds):
        cluster_allocations2[cluster_allocations == e] = i
    return counts[reorder_inds], unique_profiles[reorder_inds, :], cluster_allocations2
