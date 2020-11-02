BLEV <- function(Phi, y, weighted=F, rand_state=NULL, plot_wgts=F) {
  svd_Phi <- svd(Phi)
  U <- svd_Phi$u
  H <- U %*% t(U)
  h <- diag(H)
  prob <- h/p
  if (plot_wgts) {
    plot(prob, t="l", ylab="Sampling probability")
  }
  indices <- sample(
    x = 1:n, 
    size = m,
    replace = T,
    prob = prob
  )
  Phi_m <- Phi[indices,]
  y_m <- y[indices]
  weights <- prob[indices]
  if (weighted) {
    beta_hat <- wls_qr(Phi_m, y_m, weights)
  } else {
    beta_hat <- qr.solve(Phi_m, y_m)
  }
  y_hat <- c(Phi %*% beta_hat)
  return(
    list(
      fitted = y_hat,
      coeff = beta_hat
    )
  )
}