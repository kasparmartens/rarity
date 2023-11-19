`%not in%` <- purrr::negate(`%in%`)

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
