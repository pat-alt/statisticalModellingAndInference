wls_qr <- function(X, y, weights) {
  Phi <- diag(weights)
  X_weighted <- sqrt(qr.solve(Phi)) %*% X
  y_weighted <- sqrt(qr.solve(Phi)) %*% y
  beta <- qr.solve(X_weighted, y_weighted)
  return(beta)
}