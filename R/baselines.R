library(Rphenograph)
library(Seurat)
library(FlowSOM)

FlowSOM_clustering <- function(Y, k){
  som_map <- SOM(Y, xdim = 15, ydim = 15, rlen = 20)
  a <- BuildMST(list(map = som_map))
  metaclusters <- as.factor(metaClustering_consensus(a$map$codes, k=k))
  cl <- metaclusters[som_map$mapping[,1]]
  cl
}

Phenograph_clustering <- function(Y, k = 30){
  Rphenograph_out <- Rphenograph(Y, k=k)
  cl <- membership(Rphenograph_out[[2]])
  as.numeric(cl)
}

Phenograph_metaclustering <- function(Y, k = 10, k_final = 5){
  Rphenograph_out <- Rphenograph(Y, k=k)
  cl <- as.numeric(membership(Rphenograph_out[[2]]))
  avg_cols <- function(mat, groups) {
    data.frame(mat, group = groups) %>% 
      group_by(group) %>% 
      summarise_all(mean) %>% 
      select(-group) %>% 
      as.matrix()
  }
  cluster_means <- avg_cols(Y, cl)
  cl_meta <- cutree(hclust(dist(cluster_means)), k_final)
  cl_meta[cl]
}


Seurat_UMAP <- function(Y, seed=1){
  # use Seurat defaults for UMAP
  mat <- uwot::umap(Y, n_neighbors = 30, min_dist = 0.3, metric = "euclidean", seed=seed)
  data.frame(UMAP_1 = mat[, 1], UMAP_2 = mat[, 2])
}

