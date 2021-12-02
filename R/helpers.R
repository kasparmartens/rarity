`%not in%` <- purrr::negate(`%in%`)

plot_rarity_unique_heatmap <- function(rarity, n_cells, scale_with_number_of_cells=FALSE, subset_clusters = NULL, subset_features = NULL, cluster_columns = TRUE, ...){
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
  
  if(scale_with_number_of_cells){
    do.call("rbind", lapply(1:how_many, function(i){
      do.call("rbind", lapply(1:counts[i], function(x){
        unique_profiles[i, ]
      }))
    }))
  } else{
    mat <- do.call("rbind", lapply(1:how_many, function(i){
      unique_profiles[i, ]
    }))
  }

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

