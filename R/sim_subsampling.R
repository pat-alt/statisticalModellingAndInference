sim_subsampling <- function(Phi, y, subsample_estimator, J=1000, bias_correct=T, B=1000, ...) {
  # Compute full sample OLS estimator:
  y_hat_ols <- qr.fitted(qr.default(Phi),y)
  # Compute J predictions from subsample estimator:
  output <- rbindlist(
    lapply(
      1:J,
      function(j) {
        estimate <- subsample_estimator(Phi, y, ...) # subsample estimate
        e <- y_hat_ols - estimate$fitted # residual
        mse <- crossprod(e) # in-sample MSE
        if (bias_correct) { 
          # Bootstrap bias-correction term:
          
        }
      }
    )
  )
}