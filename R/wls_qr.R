wls_qr <- function(X, y, weights) {
  Phi <- diag(weights)
  X_weighted <- sqrt(Phi) %*% X
  y_weighted <- sqrt(Phi) %*% y
  beta <- qr.solve(X_weighted, y_weighted)
  return(beta)
}