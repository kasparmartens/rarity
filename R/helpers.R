`%not in%` <- purrr::negate(`%in%`)

plot_binary_signatures <- function(rarity, n_cells = 30, subset_clusters = NULL, subset_features = NULL, cluster_columns = TRUE, ...){
  if(is.null(subset_clusters)){
    counts <- rarity$counts
    unique_profiles <- rarity$unique_profiles
  } else{
    counts <- rarity$counts[subset_clusters]
    unique_profiles <- rarity$unique_profiles[subset_clusters, , drop=FALSE]
  }
  
  if(!is.null(subset_features)){
    unique_profiles <- unique_profiles[, subset_features, drop=FALSE]
  }
  
  how_many <- sum(counts > n_cells)
  
  mat <- do.call("rbind", lapply(1:how_many, function(i){
    unique_profiles[i, ]
  }))
  
  rownames(mat) <- sprintf("%5d cells", counts[1:how_many])

  ComplexHeatmap::Heatmap(
    mat,
    col = RColorBrewer::brewer.pal(5, "Blues"),
    row_split = subset_clusters, 
    cluster_columns = cluster_columns,
    cluster_rows = FALSE,
    show_column_dend = FALSE,
    heatmap_legend_param = list(title = ""),
    rect_gp = grid::gpar(col = "white", lwd=2),
    use_raster = TRUE, 
    show_heatmap_legend = FALSE,
    ...
  )
}

plot_heatmap_with_observations <- function(Y, rarity, n_cells = 30, subset_cells = NULL, cluster_columns = TRUE, ...){
  cluster_allocations <- rarity$cluster_allocations
  n_clusters <- max(cluster_allocations)
  tiny_clusters <- which(rarity$counts < n_cells)
  if(is.null(subset_cells)){
    subset_cells <- (cluster_allocations %not in% tiny_clusters)
  } else{
    subset_cells <- subset_cells & (cluster_allocations %not in% tiny_clusters)
  }
  
  ComplexHeatmap::Heatmap(
    Y[subset_cells, ],
    col = c(RColorBrewer::brewer.pal(9, "Reds"), "#67000D"),
    cluster_columns = cluster_columns,
    row_split = cluster_allocations[subset_cells], 
    cluster_row_slices = FALSE,
    cluster_rows = TRUE,
    use_raster = TRUE,
    heatmap_legend_param = list(title = ""),
    show_row_dend = FALSE,
    show_column_dend = FALSE, 
    ...
  )
}

boxplot_diffexp_selected_clusters <- function(mat, cl){
  df <- mat %>% 
    data.frame() %>% 
    mutate(cluster = factor(cl)) %>% 
    gather(marker, value, -cluster)
  
  df %>% 
    ggplot(aes(cluster, value, col=cluster)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ marker) +
    theme_bw() +
    theme(legend.position = "none")
}

mylog <- function(p) ifelse(p==0, 0, log(p))

entropy <- function(x){
  probs <- x / sum(x)
  -sum(probs * mylog(probs))
}


helper_metrics <- function(inferred, true){
  
  homogeneity <- clevr::homogeneity(true, inferred)
  completeness <- clevr::completeness(true, inferred)
  
  list(
    homogeneity = homogeneity,
    completeness = completeness,
    v_measure = 2 * homogeneity * completeness / (homogeneity + completeness)
  )
}

helper_metrics_conditional <- function(freq_table, selected_cell_type){
  # entropy K given C
  entropy_K_given_C <- entropy(freq_table[, colnames(freq_table) == selected_cell_type])
  entropy_K <- entropy(rep(1, nrow(freq_table))) #entropy(rowSums(freq_table))
  # entropy C given argmax K
  most_likely_k <- which.max(freq_table[, colnames(freq_table) == selected_cell_type])
  entropy_C_given_K <- entropy(freq_table[most_likely_k, ])
  entropy_C <- entropy(rep(1, ncol(freq_table))) #entropy(colSums(freq_table))
  
  homogeneity <- 1.0 - entropy_C_given_K / entropy_C
  completeness <- 1.0 - entropy_K_given_C / entropy_K
  
  list(
    homogeneity = homogeneity,
    completeness = completeness,
    v_measure = 2 * homogeneity * completeness / (homogeneity + completeness)
  )
}


Seurat_UMAP <- function(Y){
  rownames(Y) <- 1:nrow(Y)
  
  seurat_object <- Seurat::CreateSeuratObject(t(Y))
  seurat_object <- Seurat::ScaleData(seurat_object)
  
  seurat_umap <- Seurat::RunUMAP(seurat_object, features = colnames(Y), seed.use = 2, metric="euclidean")
  
  data.frame(seurat_umap@reductions$umap@cell.embeddings)
}

Phenograph_clustering <- function(Y, k = 30){
  library(Rphenograph)
  Rphenograph_out <- Rphenograph(Y, k=k)
  cl <- membership(Rphenograph_out[[2]])
  as.numeric(cl)
}

