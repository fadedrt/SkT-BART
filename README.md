# SkT-BART

This repository contains code to implement the methodology described in the paper *"Bayesian Additive Regression Trees for Responses with Heavy Tails and Skewness"*. 

The Bayesian tree structure and framework in this code are adapted from the primary functions of [CSP-BART](https://github.com/ebprado/CSP-BART/tree/main). Model comparisons are performed using the [BART](https://github.com/ebprado/CSP-BART/blob/main/cspbart/R/bart.R) module from CSP-BART and the implementation of [skewBART](https://github.com/Seungha-Um/skewBART?tab=readme-ov-file).

## Example

```r
# Load required libraries
library(skewt)   # For skewed t distribution functions
library(zeallot) # For multiple assignment (%<-%)
library(rust)    # For the update of skew parameter gamma

# Function to simulate Friedman's dataset with skew-t errors
sim_fried <- function(N, P, v, sigma, gamma) {
  X <- matrix(runif(N * P), nrow = N)  # Random uniform predictors
  mu <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5] # True mean function
  Y <- mu + sigma * rskt(N, v, gamma)  # Generate responses with skew-t errors
  return(list(X = X, Y = Y, mu = mu, v = v, gamma = gamma))
}

# -------------------------
# Simulate training and test datasets
# -------------------------
set.seed(123)
c(x, y, mu, v, gamma) %<-% sim_fried(250, 5, v = 1, sigma = 1, gamma = 1)

set.seed(123)
c(test_x, test_y, test_mu, v_t, gamma_t) %<-% sim_fried(100, 5, v = 1, sigma = 1, gamma = 1)

# Fit sktBART model
sktbart_fit <- sktbart(y, x, ntrees = 50, nburn = 2500, npost = 2500)

# Compute estimated parameters
v_hat <- mean(sktbart_fit$v)
sigma_hat <- sqrt(mean(sktbart_fit$sigma2))
gamma_hat <- sqrt(mean(sktbart_fit$gamma2))

# Generate predictions
c(f_hat_train, f_hat_test, f_hat_train_mean, f_hat_test_mean) %<-% 
  prediction(x, y, test_x, npost = 2500, trees = sktbart_fit$trees)

# Generate residual shift based on estimated parameters
if (v_hat > 1) {
  # Use the median of the sampled skewed errors for shifting to ensure numerical stability
  # particularly under heavy-tailed scenarios
  set.seed(123)
  r_shift <- median(sigma_hat * rskt(10000, v_hat, gamma_hat))
} else {
  # When degrees of freedom <= 1, expectation does not exist
  r_shift <- 0
}

# Adjust predictions
final_predictions <- f_hat_test_mean + r_shift

# -------------------------
# Evaluation metrics
# -------------------------
cat("sktBART RMSE:", round(sqrt(mean((final_predictions - test_y)^2)), 2), "\n")
cat("sktBART MAE:", round(mean(abs(final_predictions - test_y)), 2), "\n")
cat("sktBART MAPE:", round(mean(abs(final_predictions - test_y) / abs(test_y)) * 100, 2), "%\n")
cat("sktBART LPML:", calculate_lpml_skewt(y = y, sktbart_fit$bart_hat,
                                    sigma2_store = sktbart_fit$sigma2, 
                                    gamma2_store = sktbart_fit$gamma2, 
                                    v_store = sktbart_fit$v), "\n")

