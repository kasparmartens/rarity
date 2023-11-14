hamming_ball_query <- function(signature_of_interest, rarity_fit, radius, min_n_cells = 30){
    # calculate hamming distances
    hamming_dists <- rowSums(abs(t(t(rarity_fit$unique_profiles) - signature_of_interest)))

    data.frame(rarity_fit$unique_profiles) %>% 
        mutate(
            cluster_idx = 1:n(), 
            counts = rarity_fit$counts, 
            hamming_dists
        ) %>% 
        filter(hamming_dists <= radius) %>% 
        filter(counts >= min_n_cells) %>%
        # return cluster indices
        pull(cluster_idx)
}
