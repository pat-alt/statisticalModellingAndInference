# Undersampler class ----
undersampler <- function(X, y) {
  distinct_values <- unique(y)
  if(length(distinct_values)>2) {
    stop("y should have 2 distinct levels.")
  }
  grouped <- lapply(
    distinct_values, 
    function(i) {
      indices <- which(y==i)
      n <- length(indices)
      list(indices=indices, n=n)
    }
  )
  min_class <- which.min(sapply(grouped, function(i) i[[2]]))
  n_min <- grouped[[min_class]]$n
  indices_min <- grouped[[min_class]]$indices
  n_maj <- grouped[[-min_class]]$n
  indices_maj <- grouped[[-min_class]]$indices
  v <- list(
    X = X,
    y = y,
    n = n_min + n_maj,
    n_min = n_min,
    indices_min = indices_min,
    n_maj = n_maj,
    indices_maj = indices_maj
  )
  class(v) <- "undersampler"
  return(v)
}

# UNIF method: ----
UNIF.undersampler <- function(vars, weighted=F, rand_state=NULL) {
  if (!is.null(rand_state)) {
    set.seed(rand_state)
  }
  invisible(list2env(vars, envir = environment()))
  indices <- sample(indices_maj, size=n_min) # sample from majority class
  indices <- c(indices_min, indices)
  X_m <- X[indices,]
  y_m <- y[indices]
  beta_hat <- glm(y_m~X_m,family = "binomial")$coefficients
  # beta_hat <- logit_irls(X_m, y_m)$coeff
  if(!all(X[,1]==1)) {
    X <- cbind(1,X)
  }
  y_hat <- c(X %*% beta_hat)
  p_y <- exp(y_hat)/(1+exp(y_hat))
  return(
    list(
      linear_predictors = y_hat,
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
  invisible(list2env(vars, envir = environment()))
  # Sampling probabilities:
  U <- svd(X)$u
  H <- tcrossprod(U)
  h <- diag(H)
  prob <- h/ncol(X)
  # Sample:
  indices <- sample(
    x = indices_maj, 
    size = n_min,
    replace = T,
    prob = prob[indices_maj]
  )
  indices <- c(indices_min, indices)
  X_m <- X[indices,]
  y_m <- y[indices]
  weights <- prob[indices]
  # Fit:
  if (weighted) {
    beta_hat <- glm(y_m~X_m,family = "binomial", weights = weights)$coefficients
  } else {
    beta_hat <- glm(y_m~X_m,family = "binomial")$coefficients
  }
  # Predict:
  if(!all(X[,1]==1)) {
    X <- cbind(1,X)
  }
  y_hat <- c(X %*% beta_hat)
  p_y <- exp(y_hat)/(1+exp(y_hat))
  
  return(
    list(
      linear_predictors = y_hat,
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
  invisible(list2env(vars, envir = environment()))
  # Leverage scores:
  U <- svd(X)$u
  H <- tcrossprod(U)
  h <- diag(H)
  # Euclidian norms:
  predictor_len <- sapply(
    1:n,
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
  # Sample:
  indices <- sample(
    x = indices_maj, 
    size = n_min,
    replace = T,
    prob = prob[indices_maj]
  )
  indices <- c(indices_min, indices)
  X_m <- X[indices,]
  y_m <- y[indices]
  weights <- prob[indices]
  # Fit:
  if (weighted) {
    beta_hat <- glm(y_m~X_m,family = "binomial", weights = weights)$coefficients
  } else {
    beta_hat <- glm(y_m~X_m,family = "binomial")$coefficients
  }
  # Predict:
  if(!all(X[,1]==1)) {
    X <- cbind(1,X)
  }
  y_hat <- c(X %*% beta_hat)
  p_y <- exp(y_hat)/(1+exp(y_hat))
  return(
    list(
      linear_predictors = y_hat,
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
  invisible(list2env(vars, envir = environment()))
  # Leverage scores:
  U <- svd(X)$u
  H <- tcrossprod(U)
  h <- diag(H)
  # Euclidian norms:
  predictor_len <- sapply(
    1:n,
    function(i) {
      norm(as.matrix(X[i,]), type="f")
    }
  )
  # Sampling probabilities:
  prob <- sapply(predictor_len, function(i) i/sum(predictor_len))
  # Sample:
  indices <- sample(
    x = indices_maj, 
    size = n_min,
    replace = T,
    prob = prob[indices_maj]
  )
  indices <- c(indices_min, indices)
  X_m <- X[indices,]
  y_m <- y[indices]
  weights <- prob[indices]
  # Fit:
  if (weighted) {
    beta_hat <- glm(y_m~X_m,family = "binomial", weights = weights)$coefficients
  } else {
    beta_hat <- glm(y_m~X_m,family = "binomial")$coefficients
  }
  # Predict:
  if(!all(X[,1]==1)) {
    X <- cbind(1,X)
  }
  y_hat <- c(X %*% beta_hat)
  p_y <- exp(y_hat)/(1+exp(y_hat))
  return(
    list(
      linear_predictors = y_hat,
      fitted = p_y,
      coeff = beta_hat
    )
  )
}

PL = function(vars, weighted=F, rand_state=NULL) {
  UseMethod("PL")
}