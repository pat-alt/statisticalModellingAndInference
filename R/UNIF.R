UNIF <- function(Phi, y, weighted=F, rand_state=NULL) {
  if (!is.null(rand_state)) {
    set.seed(rand_state)
  }
  indices <- sample(1:n, size=m)
  Phi_m <- Phi[indices,]
  y_m <- y[indices]
  beta_hat <- qr.solve(Phi_m, y_m)
  y_hat <- c(Phi %*% beta_hat_random)
  return(
    list(
      fitted = y_hat,
      coeff = beta_hat,
    )
  )
}