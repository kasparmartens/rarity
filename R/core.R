Rarity_fit <- function(Y, n_iter = 20000L, batch_size = 2500L){
  if(!is.integer(n_iter)) n_iter <- as.integer(n_iter)
  if(!is.integer(batch_size)) batch_size <- as.integer(batch_size)
  
  my_module <- import_from_path("core", "rarity")
  out <- my_module$fit_Rarity(Y = Y, n_iter=n_iter, batch_size=batch_size)
  # postprocess output
  colnames(out$unique_profiles) <- colnames(Y)
  colnames(out$z_binary) <- colnames(Y)
  out
}
