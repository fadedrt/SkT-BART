# SkT-BART
This repository contains code to implement the methodology described in the paper “Bayesian Adve Regression Trees for Responses with Heavy Tails and Skewness”. 
The Bayesian tree structure and framework in this code are adapted from the primary functions of [CSP-BART](https://github.com/ebprado/CSP-BART/tree/main).

Model comparisons are performed using the [BART](https://github.com/ebprado/CSP-BART/blob/main/cspbart/R/bart.R) module from CSP-BART and the implementation of [skewBART](https://github.com/Seungha-Um/skewBART?tab=readme-ov-file).

## Example
```r
```r
# Load required libraries
library(skewt)   # For skewed t distribution functions
library(zeallot) # For multiple assignment (%<-%)
library(rust)    # Optional, depending on your implementation

# Function to simulate Friedman's dataset with skew-t errors
sim_fried <- function(N, P, v, sigma, gamma) {
  X <- matrix(runif(N * P), nrow = N)  # Random uniform predictors
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5] # True mean function
  Y <- mu + sigma * rskt(N, v, gamma)  # Generate responses with skew-t errors
  return(list(X = X, Y = Y, mu = mu, v = v, gamma = gamma))
}

# Function to generate predictions for train and test sets
prediction <- function(x, y, test_x, test_y, npost, trees){
  y_mean = mean(y)
  y_sd = sd(y)
  prediction_train = matrix(NA, ncol = npost, nrow = nrow(x))
  prediction_test = matrix(NA, ncol = npost, nrow = nrow(test_x))
  
  # Generate predictions for training set
  for (i in 1:npost){
    prediction_train[, i] <- get_predictions(trees[[i]], x, single_tree = FALSE)
  }
  
  # Generate predictions for test set
  for (i in 1:npost){
    prediction_test[, i] <- get_predictions(trees[[i]], test_x, single_tree = FALSE)
  }
  
  # Rescale predictions back to original scale
  f_hat_train      <- prediction_train * y_sd + y_mean
  f_hat_test       <- prediction_test * y_sd + y_mean
  f_hat_train_mean <- rowMeans(f_hat_train)
  f_hat_test_mean  <- rowMeans(f_hat_test)
  
  return(list(f_hat_train = f_hat_train, f_hat_test = f_hat_test, 
              f_hat_train_mean = f_hat_train_mean, 
              f_hat_test_mean = f_hat_test_mean))
}

# Function to calculate LPML (Log Pseudo-Marginal Likelihood)
calculate_lpml <- function(y, x, npost = 100, trees, sigma2_store, gamma2_store, v_store) {
  y_mean <- mean(y)
  y_sd <- sd(y)
  n <- length(y)
  
  # Predictive values for each MCMC sample
  prediction_trees <- matrix(NA, nrow = n, ncol = npost)
  for (j in 1:npost) {
    prediction_trees[, j] <- get_predictions(trees[[j]], x, single_tree = FALSE) * y_sd + y_mean
  }
  
  # Initialize CPO (Conditional Predictive Ordinate)
  cpo <- numeric(n)
  for (i in 1:n) {
    predictive_densities <- numeric(npost)
    for (j in 1:npost) {
      mu_ij <- prediction_trees[i, j]
      sigma_j <- sqrt(sigma2_store[j])
      gamma_j <- sqrt(gamma2_store[j])
      v_j <- v_store[j]
      
      z <- (y[i] - mu_ij) / sigma_j
      predictive_densities[j] <- dskt(z, df = v_j, gamma = gamma_j) / sigma_j
    }
    
    # Avoid numerical underflow
    predictive_densities <- pmax(predictive_densities, 1e-8)
    cpo[i] <- 1 / mean(1 / predictive_densities)  # Harmonic mean
  }
  
  # Compute LPML
  lpml <- sum(log(cpo))
  return(lpml)
}

# -------------------------
# Simulate training and test datasets
# -------------------------
set.seed(123)
c(x, y, mu, v, gamma) %<-% sim_fried(250, 5, v = 1, sigma = 1, gamma = 1)

set.seed(123)
c(test_x, test_y, test_mu, v, gamma) %<-% sim_fried(100, 5, v = 1, sigma = 1, gamma = 1)

# Fit sktBART model
sktbart_fit = sktbart(y, x, ntrees = 50, nburn = 2500, npost = 2500)

# Compute estimated parameters
v_hat = mean(sktbart_fit$v)
sigma_hat = sqrt(mean(sktbart_fit$sigma2))
gamma_hat = sqrt(mean(sktbart_fit$gamma2))

# Generate residuals based on estimated parameters
if (v_hat > 1) {
  r = sigma_hat * rskt(10000, v_hat, gamma_hat)
} else if (v_hat > 0 && v_hat <= 1) {
  # When degrees of freedom <= 1, expectation does not exist
  r = 0
}

# Generate predictions
c(f_hat_train, f_hat_test, f_hat_train_mean, f_hat_test_mean) %<-% 
  prediction(x, y, test_x, test_y, npost = 2500, trees = sktbart_fit$trees)

# -------------------------
# Evaluation metrics
# -------------------------
cat("sktBART RMSE:", round(sqrt(sum((f_hat_test_mean - test_y)^2) / length(test_y)), 2), "\n")
cat("sktBART MAE:", round(sum(abs(f_hat_test_mean - test_y)) / length(test_y), 2), "\n")
cat("sktBART MAPE:", round(mean(abs(f_hat_test_mean - test_y) / abs(test_y)) * 100, 2), "\n")
cat("sktBART LPML:", calculate_lpml(y = y, x = x, npost = 2500, 
                                    trees = sktbart_fit$trees, 
                                    sigma2_store = sktbart_fit$sigma2, 
                                    gamma2_store = sktbart_fit$gamma2, 
                                    v_store = sktbart_fit$v))

