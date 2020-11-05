# Undersampler class ----
undersampler <- function(X, y) {
  n_min <- min(table(y))
  n_maj <- max(table(y))
  if(!all(X[,1]==1)) {
    X <- cbind(1,X)
  }
  v <- list(
    X = X,
    y = y,
    n_min = n_min,
    n_maj = n_maj
  )
  class(v) <- "undersampler"
  return(v)
}

# UNIF method: ----
UNIF.undersampler <- function(vars, weighted=F, rand_state=NULL) {
  if (!is.null(rand_state)) {
    set.seed(rand_state)
  }
  X <- vars$X
  y <- vars$y
  m <- vars$n_min
  n <- m + vars$n_maj
  indices <- sample(1:n, size=m)
  X_m <- X[indices,]
  y_m <- y[indices]
  beta_hat <- logit_irls(X_m, y_m)$coeff
  y_hat <- c(X %*% beta_hat)
  p_y <- exp(y_hat)/(1+exp(y_hat))
  
  return(
    list(
      fitted = p_y,
      coeff = beta_hat
    )
  )
}

UNIF = function(vars, weighted=F, rand_state=NULL) {
  UseMethod("UNIF")
}

# BLEV method: ----
BLEV.undersampler <- function(vars, weighted=F, rand_state=NULL) {
  if (!is.null(rand_state)) {
    set.seed(rand_state)
  }
  X <- vars$X
  y <- vars$y
  m <- vars$n_min
  n <- m + vars$n_maj
  # Sampling probabilities:
  U <- svd(X)$u
  H <- tcrossprod(U)
  h <- diag(H)
  prob <- 1-(h/ncol(X))
  # Sample:
  indices <- sample(
    x = 1:n, 
    size = m,
    replace = T,
    prob = prob
  )
  X_m <- X[indices,]
  y_m <- y[indices]
  weights <- prob[indices]
  # Fit:
  if (weighted) {
    beta_hat <- logit_irls(X_m, y_m, weights)$coeff
  } else {
    beta_hat <- logit_irls(X_m, y_m)$coeff
  }
  # Predict:
  y_hat <- c(X %*% beta_hat)
  p_y <- exp(y_hat)/(1+exp(y_hat))
  
  return(
    list(
      fitted = p_y,
      coeff = beta_hat
    )
  )
}

BLEV = function(vars, weighted=F, rand_state=NULL) {
  UseMethod("BLEV")
}

# OPT method: ----
OPT.undersampler <- function(vars, weighted=F, rand_state=NULL) {
  if (!is.null(rand_state)) {
    set.seed(rand_state)
  }
  X <- vars$X
  y <- vars$y
  m <- vars$n_min
  n <- m + vars$n_maj
  # Leverage scores:
  U <- svd(X)$u
  H <- tcrossprod(U)
  h <- diag(H)
  # Euclidian norms:
  predictor_len <- sapply(
    1:nrow(X),
    function(i) {
      norm(as.matrix(X[i,]), type="f")
    }
  )
  # Optimal sampling probabilities:
  prob <- sapply(
    1:n, 
    function(i) {
      (sqrt(1-h[i]) * predictor_len[i]) / crossprod(sqrt(1-h),predictor_len)[1]
    }
  )
  prob = 1 - prob
  # Sample:
  indices <- sample(
    x = 1:n, 
    size = m,
    replace = T,
    prob = prob
  )
  X_m <- X[indices,]
  y_m <- y[indices]
  weights <- prob[indices]
  # Fit:
  if (weighted) {
    beta_hat <- logit_irls(X_m, y_m, weights)$coeff
  } else {
    beta_hat <- logit_irls(X_m, y_m)$coeff
  }
  # Predict:
  y_hat <- c(X %*% beta_hat)
  p_y <- exp(y_hat)/(1+exp(y_hat))
  
  return(
    list(
      fitted = p_y,
      coeff = beta_hat
    )
  )
}

OPT = function(vars, weighted=F, rand_state=NULL) {
  UseMethod("OPT")
}

# PL method: ----
PL.undersampler <- function(vars, weighted=F, rand_state=NULL) {
  if (!is.null(rand_state)) {
    set.seed(rand_state)
  }
  X <- vars$X
  y <- vars$y
  m <- vars$n_min
  n <- m + vars$n_maj
  # Euclidian norms:
  predictor_len <- sapply(
    1:nrow(X),
    function(i) {
      norm(as.matrix(X[i,]), type="f")
    }
  )
  # Sampling probabilities:
  prob <- sapply(predictor_len, function(i) i/sum(predictor_len))
  prob = 1 - prob
  # Sample:
  indices <- sample(
    x = 1:n, 
    size = m,
    replace = T,
    prob = prob
  )
  X_m <- X[indices,]
  y_m <- y[indices]
  weights <- prob[indices]
  # Fit:
  if (weighted) {
    beta_hat <- logit_irls(X_m, y_m, weights)$coeff
  } else {
    beta_hat <- logit_irls(X_m, y_m)$coeff
  }
  # Predict:
  y_hat <- c(X %*% beta_hat)
  p_y <- exp(y_hat)/(1+exp(y_hat))
  return(
    list(
      fitted = p_y,
      coeff = beta_hat
    )
  )
}

PL = function(vars, weighted=F, rand_state=NULL) {
  UseMethod("PL")
}