library(tidyverse)

rmatnorm <- function(mean_mat, sd_mat){
  eps <- rnorm(nrow(mean_mat) * ncol(mean_mat), mean=mean_mat, sd=sd_mat)
  matrix(eps, nrow(mean_mat), ncol(mean_mat))
}

generate_synthetic_data <- function(marker_mat, n_cells){
  assertthat::assert_that(nrow(marker_mat) == length(n_cells))
  
  Y <- do.call("rbind", lapply(1:nrow(marker_mat), function(i){
    binary_exprs <- t(replicate(n_cells[i], marker_mat[i, ]))
    epsilon <- rnorm(ncol(binary_exprs) * nrow(binary_exprs))
    mean <- 0.05 + 0.5 * binary_exprs
    sd <- 0.03 * (binary_exprs == 0) + 0.18 * (binary_exprs == 1)
    exprs <- rmatnorm(mean, sd)
    exprs[exprs < 0] <- 0
    exprs[exprs > 1] <- 0.99
    exprs
  }))
  colnames(Y) <- paste0("marker", 1:ncol(Y))
  Y
}

# binary matrix for cell types A-E and seven markers
marker_mat <- rbind(
  # cell type A
  c(1, 1, 1, 1, 0, 0, 0),
  # cell type B
  c(1, 1, 1, 0, 1, 0, 0),
  # cell type C
  c(1, 1, 1, 1, 0, 1, 0),
  # cell type D
  c(1, 1, 1, 0, 0, 0, 0),
  # cell type E
  c(1, 1, 0, 0, 0, 0, 0)
)


n_cells <- c(4000, 3000, 1000, 60, 40)
rownames(marker_mat) <- sprintf("cluster %d, # cells = %d", 1:length(n_cells), n_cells)
colnames(marker_mat) <- sprintf("gene %s", 1:ncol(marker_mat))

set.seed(1)
Y <- generate_synthetic_data(marker_mat, n_cells)
true_cl <- rep(1:nrow(marker_mat), n_cells)

write_csv(data.frame(Y), "data/synthetic1.csv")
