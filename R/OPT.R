OPT <- function(Phi, y, m, weighted=F, rand_state=NULL, plot_wgts=F, prob_only=F) {
  n <- nrow(Phi)
  # Leverage scores:
  svd_Phi <- svd(Phi)
  U <- svd_Phi$u
  H <- tcrossprod(U)
  h <- diag(H)
  # Euclidian norms:
  predictor_len <- sqrt(Phi**2 %*% rep(1,ncol(Phi)))
  # Optimal sampling probabilities:
  prob <- (sqrt(1-h) * predictor_len) / crossprod(sqrt(1-h),predictor_len)[1]
  # Plot:
  if (plot_wgts) {
    plot(prob, t="l", ylab="Sampling probability")
  }
  # Output:
  if (prob_only) {
    return(prob)
  } else {
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
        coeff = beta_hat,
        prob = prob
      )
    )
  }
}