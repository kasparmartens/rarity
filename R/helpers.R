
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


plot_UMAP_with_highlighted_clusters <- function(cl, true_cl, df_umap, shuffle_colors = FALSE, highlighted_clusters = NULL, subset_cells=1:length(true_cl)){
  tbl <- table(cl, true_cl)
  if(is.null(highlighted_clusters)){
    argmax_idx <- apply(tbl, 2, which.max)[4:5]
    highlighted_clusters <- rownames(tbl)[argmax_idx]
  }
  
  
  if(shuffle_colors){
    df_umap <- df_umap %>% 
      mutate(clusters = forcats::fct_shuffle(factor(cl)))
  } else{
    df_umap <- df_umap %>% 
      mutate(clusters = factor(cl))
  }
  
  df_umap[subset_cells, ] %>% 
    ggplot(aes(UMAP_1, UMAP_2, col=clusters)) +
    geom_point(aes(alpha=clusters %in% highlighted_clusters)) +
    scale_alpha_manual(guide="none", values = c("TRUE" = 0.75, "FALSE" = 0.1)) +
    labs(x = expression(UMAP[1]), y = expression(UMAP[2]), col="Clusters") +
    theme_classic()
}

helper_plot_rarity <- function(rarity, true_cl, df_umap, n_cells = 30){
  cl <- data.frame(cl = rarity$cluster_allocations) %>% 
    mutate(cl = case_when(
      cl == 1 ~ "A",
      cl == 2 ~ "B",
      cl == 3 ~ "C",
      TRUE ~ sprintf("cluster %s", cl)
    )) %>% 
    pull(cl)
  tiny_clusters <- which(rarity$counts < n_cells)
  selected_cells <- (rarity$cluster_allocations %not in% tiny_clusters)
  highlighted_clusters <- paste("cluster", 4:5)
  
  df_umap %>% 
    mutate(clusters = factor(cl)) %>% 
    filter(selected_cells) %>% 
    ggplot(aes(UMAP_1, UMAP_2, col=clusters)) +
    geom_point(aes(alpha = clusters %in% highlighted_clusters)) +
    scale_alpha_manual(guide="none", values = c("TRUE" = 0.75, "FALSE" = 0.1)) +
    labs(x = expression(UMAP[1]), y = expression(UMAP[2]), col="Clusters") +
    theme_classic() + 
    scale_color_viridis_d(option = "magma")
}

helper_plot_UMAP_rare_cells <- function(df_umap, true_cl){
  
  df_umap_with_labels <- df_umap %>% 
    mutate(clusters = case_when(
      true_cl == 1 ~ "cell type A",
      true_cl == 2 ~ "cell type B",
      true_cl == 3 ~ "cell type C",
      true_cl == 4 ~ "cell type D \n(rare)",
      true_cl == 5 ~ "cell type E (rare)"
    ))
  
  df_umap_with_labels_summarised <- df_umap_with_labels %>% 
    group_by(clusters) %>% 
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) %>% 
    mutate(UMAP_2 = case_when(
      clusters == "cell type B" ~ UMAP_2 + 1,
      clusters == "cell type D \n(rare)" ~ UMAP_2 - 2.5,
      clusters == "cell type E (rare)" ~ UMAP_2 + 1.7,
      TRUE ~ UMAP_2
    ))
  
  df_umap_with_labels %>% 
    ggplot(aes(UMAP_1, UMAP_2, col=clusters, alpha = clusters)) +
    geom_point() +
    geom_label(aes(label=clusters), data = df_umap_with_labels_summarised, alpha=0.9, size=5, show.legend = FALSE) +
    scale_color_viridis_d(direction = -1, end=0.9, option="magma") +
    scale_alpha_manual(values = c("cell type A" = 0.4,
                                  "cell type B" = 0.3,
                                  "cell type C" = 0.3,
                                  "cell type D \n(rare)" = 0.5,
                                  "cell type E (rare)" = 0.75)) +
    labs(x = expression(UMAP[1]), y = expression(UMAP[2]), col="Clusters", alpha="Clusters") +
    theme_classic() +
    theme(legend.position = "none")
}
