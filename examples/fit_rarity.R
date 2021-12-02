library(tidyverse)
library(reticulate)
use_python("/usr/local/bin/python3")

source("R/core.R")
source("R/helpers.R")

# read in the synthetic data
df <- read_csv("data/synthetic1.csv")
Y <- as.matrix(df)

# train the model
m <- Rarity_fit(Y, n_iter = 20000L, batch_size = 512)

# plot unique binary profiles
plot_rarity_unique_heatmap(m, n_cells = 30)

