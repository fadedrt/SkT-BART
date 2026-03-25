#' MCMC Update for Skew-t Error Distribution Parameters
#' 
#' @param x Residuals from the BART model
#' @param gamma2 Current skewness parameter (square)
#' @param v Current degrees of freedom
#' @param d Hyperparameter for v update (default 0.01)
#' @param alpha_gamma Hyperparameter for gamma2 prior
#' @param beta_gamma Hyperparameter for gamma2 prior
#' @param lambda Latent weight variables (Data Augmentation)
#' @param sigma_alpha Hyperparameter for sigma2 Inverse-Gamma prior (shape)
#' @param sigma_beta Hyperparameter for sigma2 Inverse-Gamma prior (rate)
#' 
#' @export
skewt_parameter_update <- function(x, gamma2, v, d = 0.1, 
                                   alpha_gamma = 2, beta_gamma = 1, 
                                   lambda, sigma_alpha, sigma_beta) {
  n <- length(x)
  
  # Indicator logic: Adjust scale based on the sign of residuals
  indicator <- ifelse(x >= 0, 1 / gamma2, gamma2)
  
  # Calculate sufficient statistic S for sigma2 update
  S <- sum(lambda * x^2 * indicator)
  
  # 1. Update Precision (sigma2)
  sigma2 <- update_sigma2_prior(S = S, n = n, 
                                sigma_alpha = sigma_alpha, 
                                sigma_beta = sigma_beta)
  
  # 2. Update Degrees of Freedom (v)
  v <- update_v_skewt(n = n, d = d, lambda1 = lambda)
  
  # 3. Update Latent Weights (lambda)
  lambda <- update_lambda_skewt(n = n, residuals = x, v = v, 
                                gamma2 = gamma2, sigma2 = sigma2)
  
  # 4. Update Skewness (gamma2)
  gamma2 <- update_gamma2_skewt(n = n, lambda = lambda, residuals = x, 
                                sigma2 = sigma2, a = alpha_gamma, b = beta_gamma)
  
  return(list(
    v = v,
    gamma2 = gamma2,
    lambda = lambda,
    sigma2 = sigma2
  ))
}

#' Update Error Variance (sigma2) with Inverse-Gamma Prior
update_sigma2_prior <- function(S, n, sigma_alpha, sigma_beta) {
  # Posterior follows Inverse-Gamma(n/2 + alpha, S/2 + beta)
  return(1 / rgamma(1, shape = n / 2 + sigma_alpha, rate = S / 2 + sigma_beta))
}

#' Update Skewness Parameter gamma2 using Ratio-of-Uniforms Sampling

update_gamma2_skewt <- function(n, lambda, residuals, sigma2, a, b) {
  log_density <- function(gamma2) {
    gamma2 <- abs(gamma2)
    indicate <- ifelse(residuals >= 0, 1 / gamma2, gamma2)
    S <- sum(lambda * residuals^2 * indicate / sigma2)
    log_prior <- (n / 2 + a - 1) * log(gamma2) - n * log(gamma2 + 1)
    log_likelihood <- -S / 2 - b * gamma2
    return(log_prior + log_likelihood)
  }
  
  lambda1 <- find_lambda_one_d(log_density)
  
  r <- ru(
    logf = log_density,
    d = 1,
    n = 1,
    trans = "BC",
    lambda = lambda1,
    lower = 0
  )
  
  return(mean(r$sim_vals)) 
}


#' Update Latent Weights (lambda)
update_lambda_skewt <- function(n, residuals, v, gamma2, sigma2) {
  # Layered representation of t-distribution via scale mixture of Normals
  rates <- (v / 2) + (residuals^2 / (2 * sigma2)) * ifelse(residuals >= 0, 1/gamma2, gamma2)
  return(rgamma(n, shape = (v + 1) / 2, rate = rates))
}


update_v_skewt <- function(n, d, lambda1){
  # Use Laplace-based mode-finding for the degrees of freedom (v).
  # The mode x_star serves as the proposal center for Rejection Sampling.
  eta = sum(lambda1 - log(lambda1))/2+d
  T=n
  f <- function(x) {
    (T / 2) * (log(x / 2) + 1 - digamma(x / 2)) + (1 / x) - eta
  }
  f_lower = f(1e-10)
  f_upper = f(100)
  # Numerical stability check: if the log-posterior is too flat or non-convex
  # at the boundaries, we use a robust empirical prior mode (v = 5).
  # This ensures MCMC stability in early iterations or with sparse data.
  if(is.nan(f_lower) || is.nan(f_upper) || f_lower*f_upper>0){
    x_star = 5
  }else{
    x_star = uniroot(f,c(1e-10,100))$root
  }
  alpha <- 1 / x_star  
  target_density <- function(x) {
    (x/2)^(T*x/2) * exp(-T * lgamma(x/2)) * exp(-eta * x)
  }
  
  proposal_density <- function(x) {
    alpha * exp(-alpha * x)
  }
  
  Q <- function(x) {
    exp((T*x/2) * log(x/2) - T * lgamma(x/2)- (T*x_star/2) * log(x_star/2) + T * lgamma(x_star/2)+ (x_star-x ) * eta-1+x/x_star)
  }
  
  ars_sample <- function(n) {
    samples <- numeric(n)
    accepted <- 0
    
    while (accepted < n) {
      x <- rexp(1, rate = alpha)
      
      prob <- Q(x)
      
      if (runif(1) < prob) {
        accepted <- accepted + 1
        samples[accepted] <- x
      }
    }
    
    return(samples)
  }
  
  return(ars_sample(1))
}


update_s_  <- function(var_count, p, alpha_s) {
  shape <- alpha_s / p + var_count
  temp  <- rgamma(length(shape), shape, rate = 1)
  temp / sum(temp)
}

                               
