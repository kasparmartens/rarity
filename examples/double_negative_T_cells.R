library(tidyverse)
library(ComplexHeatmap)

source("R/core.R")
source("R/helpers.R")

# first download the source file from zenodo https://doi.org/10.5281/zenodo.6029530
data <- read_csv("data/colon_mucosa_dataset.csv")

# data matrix
Y <- data %>%
  select(CD45, CD3, CD4, CD8, CD45RA, CD45RO, FoxP3, PD1, Ki67) %>%
  as.matrix()

# UMAP
df_umap <- Seurat_UMAP(Y)

# # fit rarity model (we can skip this step and directly download the output file)
# library(reticulate)
# use_python("/usr/local/bin/python3")
# m <- Rarity_fit(Y, n_iter = 20000L, batch_size = 1000L)
# # visualise all T-cells
# is_Tcell <- (m$z_binary[, "CD45"] > 0) & (m$z_binary[, "CD3"] > 0)
# is_Tcell_cluster <- (m$unique_profiles[, "CD45"] > 0) & (m$unique_profiles[, "CD3"] > 0)
# double_negatives <- (m$z_binary[, "CD45"] > 0) & (m$z_binary[, "CD3"] > 0) & (m$z_binary[, "CD4"] == 0) & (m$z_binary[, "CD8"] == 0)
# double_negatives_cluster <- (m$unique_profiles[, "CD45"] > 0) & (m$unique_profiles[, "CD3"] > 0) & (m$unique_profiles[, "CD4"] == 0) & (m$unique_profiles[, "CD8"] == 0)
# df_output = data.frame(
#   cell_idx = 1:length(is_Tcell), 
#   cluster_idx = m$cluster_allocations,
#   is_Tcell = is_Tcell, 
#   is_double_negative = double_negatives
# )

# download the saved output file from zenodo https://doi.org/10.5281/zenodo.6029530
df_output <- read_csv("output/rarity_results.csv") %>% 
  mutate(cluster_idx = ifelse(cluster_idx == "other", 151, as.numeric(cluster_idx)))
is_Tcell <- df_output$is_Tcell

# subset corresponding to T-cells
df_output_Tcells <- df_output %>% 
  filter(is_Tcell)

### UMAP on T-cells

df_umap <- Seurat_UMAP(Y[is_Tcell, ])

cl <- Phenograph_clustering(Y[is_Tcell, ], k=30)

df_umap %>% 
  bind_cols(select(data[is_Tcell, ], CD4, CD8, FoxP3)) %>% 
  ggplot(aes(UMAP_1, UMAP_2, col=CD4)) +
  geom_point(alpha=0.1) +
  scale_color_viridis_c() +
  labs(x = expression(UMAP[1]), y = expression(UMAP[2])) +
  theme_classic() -> p1

df_umap %>% 
  bind_cols(select(data[is_Tcell, ], CD4, CD8, FoxP3)) %>% 
  ggplot(aes(UMAP_1, UMAP_2, col=CD8)) +
  geom_point(alpha=0.1) +
  scale_color_viridis_c() +
  labs(x = expression(UMAP[1]), y = expression(UMAP[2])) +
  theme_classic() -> p2

df_umap %>% 
  bind_cols(select(data[is_Tcell, ], CD4, CD8, FoxP3)) %>% 
  ggplot(aes(UMAP_1, UMAP_2, col=FoxP3)) +
  geom_point(alpha=0.1) +
  scale_color_viridis_c() +
  labs(x = expression(UMAP[1]), y = expression(UMAP[2])) +
  theme_classic() -> p3

df_umap %>% 
  mutate(cluster = factor(
    ifelse(
      df_output_Tcells$cluster_idx %in% c(80, 94, 127, 128), df_output_Tcells$cluster_idx, "Other T cell"
    )
  )) %>% 
  arrange(desc(cluster)) %>% 
  mutate(cluster = forcats::fct_relevel(factor(cluster), c(80, 94, 127, 128, "Other T cell"))) %>% 
  ggplot(aes(UMAP_1, UMAP_2, col=cluster)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values = c("Other T cell" = "grey90", "80" = "#a50f15", "94" = "#cb181d", "127" = "#fb6a4a", "128" = "#fcbba1")) +
  theme_classic() +
  labs(x = expression(UMAP[1]), y = expression(UMAP[2])) +
  labs(col = "Double-negative\nT cell clusters") -> p4


p <- p1 + p2 + p3 + p4 + patchwork::plot_layout(nrow=2)
ggsave("figures/fig_double_negative_UMAP.png", p, width=7.4, height=5.4)


#### Heatmap

plot_rarity_heatmap_without_sample_labels <- function(Y, rarity, sample_IDs, subset_cells = 1:n_cells, n_cells = 10){
  # cluster_allocations <- rarity$cluster_allocations
  # n_clusters <- max(cluster_allocations)
  # tiny_clusters <- which(rarity$counts < n_cells)
  # selected_cells <- (cluster_allocations %not in% tiny_clusters) & subset_cells
  cluster_allocations <- asdf
  counts <- df_output %>% group_by(cluster_idx) %>% count()
  n_clusters <- max(df_cl)
  
  ComplexHeatmap::Heatmap(
    Y[selected_cells, ],
    col = c(RColorBrewer::brewer.pal(9, "Reds"), "#67000D"),
    cluster_columns = TRUE,
    row_split = cluster_allocations[selected_cells], 
    cluster_row_slices = FALSE,
    cluster_rows = TRUE,
    use_raster = TRUE,
    heatmap_legend_param = list(title = ""),
    show_column_dend = FALSE
  )
}

# produce the heatmap figure
png("figures/fig_double_negative_heatmap.png", width=115, height=180, units="mm", res = 300)

selected_cells <- df_output$cluster_idx %in% c(80, 94, 127, 128)
ComplexHeatmap::Heatmap(
  Y[selected_cells, ],
  col = c(RColorBrewer::brewer.pal(9, "Reds"), "#67000D"),
  cluster_columns = TRUE,
  row_split = df_output$cluster_idx[selected_cells], 
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  use_raster = TRUE,
  heatmap_legend_param = list(title = ""),
  show_column_dend = FALSE
)

dev.off()  


### alluvial

library(ggalluvial)

# phenograph clustering
cl_phenograph <- Phenograph_clustering(Y[is_Tcell, ])

df_alluvial <- data.frame(
  rarity = df_output$cluster_idx[is_Tcell],
  phenograph = cl_phenograph
) %>% 
  group_by(rarity, phenograph) %>% 
  summarise(freq = n()) %>% 
  mutate(
    is_double_negative = factor(
      case_when(
        rarity %in% c(80, 94, 127, 128) ~ as.character(rarity),
        TRUE ~ "Other T cell"
      )
    ),
    double_negative_cluster = forcats::fct_relevel(is_double_negative, c("Other T cell", "127", "128", "80", "94")),
    NULL
  )

df_alluvial %>% 
  filter(rarity %in% 1:130) %>%
  arrange(desc(rarity)) %>% 
  mutate(double_negative_cluster = forcats::fct_relevel(factor(double_negative_cluster), c(80, 94, 127, 128, "Other T cell"))) %>% 
  ggplot(aes(y=freq, axis1=rarity, axis2=phenograph)) +
  geom_flow(aes(fill=double_negative_cluster, col=double_negative_cluster), stat = "flow", alpha=0.75) +
  geom_stratum(aes()) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), min.y=30) +
  scale_fill_manual(values = c("Other T cell" = "grey92", "80" = "#a50f15", "94" = "#cb181d", "127" = "#fb6a4a", "128" = "#fcbba1")) +
  scale_color_manual(values = c("Other T cell" = "grey85", "80" = "#a50f15", "94" = "#cb181d", "127" = "#fb6a4a", "128" = "#fcbba1")) +
  scale_x_continuous(breaks = 1:2, labels = c("Rarity", "Phenograph")) +
  scale_y_continuous(breaks=NULL) +
  theme_void() +
  labs(y = "", col = "Double-negative\nT cell cluster", fill="Double-negative\nT cell cluster")

ggsave("figures/fig_double_negative_alluvial.png", width=5.6, height=6.3)
