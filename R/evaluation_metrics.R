# ==============================================================================
# Script: evaluation_metrics.R
# Purpose: Utility functions for generating predictions and evaluating model fit 
#          (RMSE, MAE, MAPE, LPML) for the SkT-BART model.
# ==============================================================================

#' Generate In-Sample and Out-of-Sample Predictions
#' 
#' @param x Training predictor matrix.
#' @param y Training response vector.
#' @param test_x Testing predictor matrix.
#' @param npost Number of posterior samples.
#' @param trees List of posterior tree ensembles.
#' @return A list containing full posterior predictions and their means.
#' @export
prediction <- function(x, y, test_x, npost, trees) {
  y_mean <- mean(y)
  y_sd   <- sd(y)
  
  prediction_train <- matrix(NA, ncol = npost, nrow = nrow(x))
  prediction_test  <- matrix(NA, ncol = npost, nrow = nrow(test_x))
  
  # Generate predictions using the tree ensembles
  for (i in 1:npost) {
    prediction_train[, i] <- get_predictions(trees[[i]], x, single_tree = FALSE)
    prediction_test[, i]  <- get_predictions(trees[[i]], test_x, single_tree = FALSE)
  }
  
  # Rescale predictions back to the original response scale
  f_hat_train <- prediction_train * y_sd + y_mean
  f_hat_test  <- prediction_test * y_sd + y_mean
  
  return(list(
    f_hat_train      = f_hat_train, 
    f_hat_test       = f_hat_test, 
    f_hat_train_mean = rowMeans(f_hat_train), 
    f_hat_test_mean  = rowMeans(f_hat_test)
  ))
}

#' Calculate Log Pseudo-Marginal Likelihood (LPML) for Skew-t Errors
#' 
#' @description Uses the Log-Sum-Exp trick for numerical stability when 
#' computing the harmonic mean of predictive densities.
#' 
#' @param y Original response vector.
#' @param f_hat_train Matrix of posterior in-sample predictions (npost x n).
#' @param sigma2_store Posterior samples of the variance parameter.
#' @param gamma2_store Posterior samples of the skewness parameter (squared).
#' @param v_store Posterior samples of the degrees of freedom.
#' @return The scalar LPML value.
#' @importFrom skewt dskt
#' @export
calculate_lpml_skewt <- function(y, f_hat_train, sigma2_store, gamma2_store, v_store) {
  npost <- nrow(f_hat_train)
  n <- length(y)
  log_cpo <- numeric(n)
  
  for (i in 1:n) {
    log_pred <- numeric(npost)
    for (s in 1:npost) {
      mu_resid <- y[i] - f_hat_train[s, i]
      sigma    <- sqrt(sigma2_store[s])
      gamma    <- sqrt(gamma2_store[s])
      v        <- v_store[s]
      
      # Compute log-density of the Skew-t distribution
      log_pred[s] <- log(dskt(mu_resid / sigma, v, gamma) / sigma)
    }
    
    # Harmonic mean calculation via the Log-Sum-Exp trick to prevent underflow
    u <- -log_pred
    max_u <- max(u)
    log_cpo[i] <- -(max_u + log(mean(exp(u - max_u))))
  }
  
  lpml <- sum(log_cpo)
  return(lpml)
}
