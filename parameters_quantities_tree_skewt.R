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
skewt_parameter_update <- function(x, gamma2, v, d = 0.01, 
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

#' Update Degrees of Freedom (v) via Accept-Reject Sampling
update_v_skewt <- function(n, d, lambda1) {
  eta <- sum(lambda1 - log(lambda1)) / 2 + d
  n_obs <- n
  
  # Solve for the optimal proposal parameter (x_star)
  target_root_fn <- function(x) {
    (n_obs / 2) * (log(x / 2) + 1 - digamma(x / 2)) + (1 / x) - eta
  }
  
  f_lower <- target_root_fn(1e-10)
  f_upper <- target_root_fn(100)
  
  if (is.nan(f_lower) || is.nan(f_upper) || f_lower * f_upper > 0) {
    x_star <- 5
  } else {
    x_star <- uniroot(target_root_fn, c(1e-10, 100))$root
  }
  
  alpha <- 1 / x_star
  
  calc_acceptance_prob <- function(x) {
    val <- (n_obs * x / 2) * log(x / 2) - 
           n_obs * lgamma(x / 2) - 
           (n_obs * x_star / 2) * log(x_star / 2) + 
           n_obs * lgamma(x_star / 2) + 
           (x_star - x) * eta - 1 + x / x_star
    return(exp(val))
  }
  
  ars_sample <- function(n_samples) {
    samples <- numeric(n_samples)
    accepted <- 0
    while (accepted < n_samples) {
      x_cand <- rexp(1, rate = alpha)
      prob <- calc_acceptance_prob(x_cand)
      if (runif(1) < prob) {
        accepted <- accepted + 1
        samples[accepted] <- x_cand
      }
    }
    return(samples)
  }
  
  return(ars_sample(1))
}

#' Laplace Approximation for Non-Exponential Family Marginal Likelihood
#' @export
tree_full_skewt_laplace <- function(tree, R, lambda, sigma2, sigma2_mu, gamma2, 
                                    common_vars, aux_factor_var) {
  
  # Step 1: Identify terminal nodes
  terminal_nodes <- as.numeric(which(tree$tree_matrix[, 'terminal'] == 1))
  
  # Step 2: Identify nodes violating factor constraints
  if (nrow(tree$tree_matrix) != 1) {
    terminal_ancestors <- get_ancestors(tree)
    ancestor_table <- table(terminal_ancestors[, 1], terminal_ancestors[, 2])
    
    terminals_invalid <- NULL
    for (k in 1:nrow(ancestor_table)) {
      ancestor_vars <- names(ancestor_table[k, ])[ancestor_table[k, ] != 0]
      is_invalid <- any(sapply(aux_factor_var, function(x) all(ancestor_vars %in% x)))
      if (is_invalid) {
        terminals_invalid[k] <- rownames(ancestor_table)[k]
      }
    }
    terminals_invalid <- as.numeric(terminals_invalid[!is.na(terminals_invalid)])
  } else {
    terminals_invalid <- 0
  }
  
  # Step 3: Handle ID mapping and prior variance alignment
  sigma2_mu_named <- rep(sigma2_mu, length(terminal_nodes))
  names(sigma2_mu_named) <- as.character(terminal_nodes)
  if (length(terminals_invalid) > 0 && any(terminals_invalid != 0)) {
    sigma2_mu_named[as.character(terminals_invalid)] <- 0
  }
  
  leaf_ids <- as.character(tree$node_indices)
  
  # Step 4: IRLS to find the local posterior mode (Laplace Approximation)
  mu_hat <- rep(0, length(R))  
  tol <- 1e-4
  max_iter <- 10  
  
  for (iter in 1:max_iter) {
    mu_old <- mu_hat
    C_i_iter <- ifelse(R > mu_hat, 1 / gamma2, gamma2)
    W_i_iter <- C_i_iter * lambda
    
    S_W_leaf  <- tapply(W_i_iter, leaf_ids, sum)
    S_WR_leaf <- tapply(W_i_iter * R, leaf_ids, sum)
    
    sigma2_mu_j <- sigma2_mu_named[names(S_W_leaf)]
    denom_mu <- S_W_leaf + (sigma2 / sigma2_mu_j)
    denom_mu[sigma2_mu_j == 0] <- Inf 
    
    mu_mode_leaf <- S_WR_leaf / denom_mu
    mu_hat <- as.numeric(mu_mode_leaf[leaf_ids])
    
    if (max(abs(mu_hat - mu_old), na.rm = TRUE) < tol) break
  }
  
  # Step 5: Finalize statistics at converged mode
  C_i <- ifelse(R > mu_hat, 1 / gamma2, gamma2)
  W_i <- C_i * lambda
  
  S_W   <- tapply(W_i, leaf_ids, sum)
  S_WR  <- tapply(W_i * R, leaf_ids, sum)
  S_WR2 <- tapply(W_i * R^2, leaf_ids, sum)
  
  sigma_mu_j_final <- sigma2_mu_named[names(S_W)]
  
  # Step 6: Analytical log-marginal likelihood
  denom <- S_W * sigma_mu_j_final + sigma2
  term1 <- log(sigma2) - log(denom) 
  term2 <- (sigma_mu_j_final * S_WR^2) / (sigma2 * denom)
  term3 <- - S_WR2 / sigma2  
  
  return(0.5 * sum(term1 + term2 + term3))
}

#' Direct Laplace Approximate Sampling for Leaf Node Parameters
#' 
#' @description Updates the 'mu' parameters by sampling from a Gaussian 
#' approximation centered at the posterior mode.
#' @export
simulate_mu_skew_laplace <- function(tree, R, lambda1, gamma2, sigma2, sigma2_mu, 
                                     common_vars, aux_factor_var) {
  
  # 1. Identify terminal nodes and constraints
  which_terminal <- as.numeric(which(tree$tree_matrix[, 'terminal'] == 1))
  
  if(nrow(tree$tree_matrix) != 1) {
    terminal_ancestors <- get_ancestors(tree)
    aux_table <- table(terminal_ancestors[, 1], terminal_ancestors[, 2])
    which_terminal_no_double_split <- NULL
    for (k in 1:nrow(aux_table)) {
      split_var_ancestors <- names(aux_table[k, ])[aux_table[k, ] != 0]
      invalid_ancestors <- any(unlist(lapply(aux_factor_var, function(x) all(split_var_ancestors %in% x))))
      if (invalid_ancestors) which_terminal_no_double_split[k] <- rownames(aux_table)[k]
    }
    which_terminal_no_double_split <- as.numeric(which_terminal_no_double_split[!is.na(which_terminal_no_double_split)]) 
  } else {
    which_terminal_no_double_split <- 0
  }
  
  sigma2_mu_aux <- rep(sigma2_mu, length(which_terminal))
  sigma2_mu_aux[which(which_terminal %in% which_terminal_no_double_split)] <- 0
  
  node_sizes <- as.numeric(tree$tree_matrix[which_terminal, 'node_size'])
  mu_values  <- numeric(length(node_sizes))
  unique_leaf_indices <- sort(unique(tree$node_indices))
  gamma2_safe <- max(gamma2, 1e-10) 
  
  # 2. Loop through each leaf to sample mu via IRLS mode-finding
  for(i in 1:length(node_sizes)) {
    
    if (sigma2_mu_aux[i] == 0) {
      mu_values[i] <- 0
      next
    }
    
    ind <- which(tree$node_indices == unique_leaf_indices[i])
    R_leaf <- R[ind]
    lambda_leaf <- lambda1[ind]
    
    valid_mask <- is.finite(R_leaf) & is.finite(lambda_leaf)
    R_leaf <- R_leaf[valid_mask]; lambda_leaf <- lambda_leaf[valid_mask]
    
    if (length(R_leaf) == 0) {
      mu_values[i] <- rnorm(1, 0, sqrt(sigma2_mu))
      next
    }
    
    # IRLS Mode Finding
    sum_lam <- max(sum(lambda_leaf), 1e-10)
    mu_hat <- sum(lambda_leaf * R_leaf) / (sum_lam + sigma2 / sigma2_mu)
    
    tolerance <- 1e-4
    for(iter in 1:10) { 
      mu_old <- mu_hat
      W_i <- ifelse(R_leaf > mu_hat, 1/gamma2_safe, gamma2_safe) * lambda_leaf
      mu_hat <- sum(W_i * R_leaf) / (sum(W_i) + sigma2 / sigma2_mu)
      if(abs(mu_hat - mu_old) < tolerance) break
    }
    
    if (!is.finite(mu_hat)) mu_hat <- median(R_leaf)
    
    # Laplace Sampling
    W_final <- ifelse(R_leaf > mu_hat, 1/gamma2_safe, gamma2_safe) * lambda_leaf
    laplace_var <- sigma2 / (sum(W_final) + sigma2 / sigma2_mu)
    laplace_sd <- sqrt(max(laplace_var, 1e-14))
    
    mu_star <- rnorm(1, mean = mu_hat, sd = laplace_sd)
    mu_values[i] <- ifelse(is.finite(mu_star), mu_star, mu_hat)
  }
  
  # 3. Map back to tree matrix
  tree$tree_matrix[, 'mu'] <- NA
  tree$tree_matrix[which_terminal, 'mu'] <- mu_values
  tree$tree_matrix[which_terminal_no_double_split, 'mu'] <- 0 
  
  return(tree)
}



