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
  v <- update_v_skewt(n = n, d = d, lambda = lambda)
  
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
  # Posterior of sigma2 follows Inverse-Gamma(n/2 + alpha, S/2 + beta)
  # Sampling from Gamma and taking the reciprocal
  return(1 / rgamma(1, shape = n / 2 + sigma_alpha, rate = S / 2 + sigma_beta))
}

#' Update Skewness Parameter gamma2 using Ratio-of-Uniforms Sampling
update_gamma2_skewt <- function(n, lambda, residuals, sigma2, a, b) {
  
  # Log-posterior density function for gamma2
  log_density <- function(g2) {
    g2 <- abs(g2)
    if(g2 <= 0) return(-1e10)
    
    indicator <- ifelse(residuals >= 0, 1 / g2, g2)
    S_scaled <- sum(lambda * residuals^2 * indicator / sigma2)
    
    # Log-Prior + Log-Likelihood
    log_prior <- (n / 2 + a - 1) * log(g2) - n * log(g2 + 1)
    log_likelihood <- -S_scaled / 2 - b * g2
    
    return(log_prior + log_likelihood)
  }
  
  # Estimate lambda parameter for ru algorithm
  l_param <- find_lambda_one_d(log_density)
  
  # Perform Ratio-of-Uniforms sampling (Single draw for MCMC validity)
  r_sample <- ru(
    logf = log_density,
    d = 1,
    n = 1, 
    trans = "BC",
    lambda = l_param,
    lower = 0
  )
  
  return(as.numeric(r_sample$sim_vals))
}

#' Update Latent Weights (lambda) - Layered Representation of t-distribution
update_lambda_skewt <- function(n, residuals, v, gamma2, sigma2) {
  # Vectorized calculation for efficiency
  # Based on the conditional Gamma distribution of the Skew-t scale mixture
  rates <- (v / 2) + (residuals^2 / (2 * sigma2)) * ifelse(residuals >= 0, 1/gamma2, gamma2)
  
  return(rgamma(n, shape = (v + 1) / 2, rate = rates))
}

update_v_skewt <- function(n, d, lambda1) {
  # Pre-compute constants
  eta <- sum(lambda1 - log(lambda1)) / 2 + d
  n_obs <- n
  
  # Define the function to solve for x_star (optimal proposal parameter)
  target_root_fn <- function(x) {
    (n_obs / 2) * (log(x / 2) + 1 - digamma(x / 2)) + (1 / x) - eta
  }
  
  # Solve for x_star using uniroot within a safe range
  f_lower <- target_root_fn(1e-10)
  f_upper <- target_root_fn(100)
  
  if (is.nan(f_lower) || is.nan(f_upper) || f_lower * f_upper > 0) {
    x_star <- 5
  } else {
    x_star <- uniroot(target_root_fn, c(1e-10, 100))$root
  }
  
  # Set the proposal distribution parameter
  alpha <- 1 / x_star
  
  # Define the log-acceptance ratio function Q(x)
  # Derived from target_density / proposal_density normalized by the value at x_star
  calc_acceptance_prob <- function(x) {
    val <- (n_obs * x / 2) * log(x / 2) - 
           n_obs * lgamma(x / 2) - 
           (n_obs * x_star / 2) * log(x_star / 2) + 
           n_obs * lgamma(x_star / 2) + 
           (x_star - x) * eta - 1 + x / x_star
    return(exp(val))
  }
  
  # Accept-Reject Sampling implementation
  ars_sample <- function(n_samples) {
    samples <- numeric(n_samples)
    accepted <- 0
    
    while (accepted < n_samples) {
      # Draw from exponential proposal distribution
      x_cand <- rexp(1, rate = alpha)
      
      # Calculate acceptance probability
      prob <- calc_acceptance_prob(x_cand)
      
      # Acceptance decision
      if (runif(1) < prob) {
        accepted <- accepted + 1
        samples[accepted] <- x_cand
      }
    }
    
    return(samples)
  }
  
  # Return a single sample for the Gibbs step
  return(ars_sample(1))
}




#' Laplace Approximation for Non-Exponential Family Marginal Likelihood
tree_full_skewt_laplace <- function(tree, R, lambda, sigma2, sigma2_mu, gamma2, 
                                    common_vars, aux_factor_var) {
  
  # -----------------------------------------------------------------------
  # Step 1: Identify all terminal nodes (leaves)
  # -----------------------------------------------------------------------
  terminal_nodes <- as.numeric(which(tree$tree_matrix[, 'terminal'] == 1))
  
  # -----------------------------------------------------------------------
  # Step 2: Identify invalid nodes based on factor constraints
  # -----------------------------------------------------------------------
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
  
  # -----------------------------------------------------------------------
  # Step 3: Create named sigma2_mu vector to handle ID mapping
  # -----------------------------------------------------------------------
  sigma2_mu_named <- rep(sigma2_mu, length(terminal_nodes))
  names(sigma2_mu_named) <- as.character(terminal_nodes)
  
  if (length(terminals_invalid) > 0 && any(terminals_invalid != 0)) {
    sigma2_mu_named[as.character(terminals_invalid)] <- 0
  }
  
  leaf_ids <- as.character(tree$node_indices)
  
  # -----------------------------------------------------------------------
  # Step 4: Laplace Approximation - IRLS to find the local posterior mode
  # -----------------------------------------------------------------------
  # Initialize mu predictions for each node to 0
  mu_hat <- rep(0, length(R))  
  
  # Iteratively find the mode with adaptive convergence check
  tol <- 1e-4
  max_iter <- 10  # Safe upper bound for iterations
  
  for (iter in 1:max_iter) {
    mu_old <- mu_hat # Store previous state
    
    # Calculate asymmetric weights based on the current mu_hat
    C_i_iter <- ifelse(R > mu_hat, 1 / gamma2, gamma2)
    W_i_iter <- C_i_iter * lambda
    
    # Compute leaf-wise sufficient statistics (weighted counts and residuals)
    S_W_leaf  <- tapply(W_i_iter, leaf_ids, sum)
    S_WR_leaf <- tapply(W_i_iter * R, leaf_ids, sum)
    
    # Align prior variances
    sigma2_mu_j <- sigma2_mu_named[names(S_W_leaf)]
    
    # mu_mode = S_WR / (S_W + sigma2 / sigma2_mu)
    denom_mu <- S_W_leaf + (sigma2 / sigma2_mu_j)
    denom_mu[sigma2_mu_j == 0] <- Inf # Handle invalid or stump nodes
    
    mu_mode_leaf <- S_WR_leaf / denom_mu
    
    # Map updated modes back to the full sample vector
    mu_hat <- as.numeric(mu_mode_leaf[leaf_ids])
    
    # Convergence check: break if the maximum update is below tolerance
    if (max(abs(mu_hat - mu_old), na.rm = TRUE) < tol) {
      break
    }
  }
  
  # -----------------------------------------------------------------------
  # Step 5: Finalize sufficient statistics based on converged mode
  # -----------------------------------------------------------------------
  C_i <- ifelse(R > mu_hat, 1 / gamma2, gamma2)
  W_i <- C_i * lambda
  
  S_W   <- tapply(W_i, leaf_ids, sum)       # Corresponds to M_j (including C_i)
  S_WR  <- tapply(W_i * R, leaf_ids, sum)   # Corresponds to S_j1
  S_WR2 <- tapply(W_i * R^2, leaf_ids, sum) # Corresponds to S_j2
  
  sigma_mu_j_final <- sigma2_mu_named[names(S_W)]
  
  # -----------------------------------------------------------------------
  # Step 6: Calculate log-marginal likelihood (Laplace integral solution)
  # -----------------------------------------------------------------------
  denom <- S_W * sigma_mu_j_final + sigma2
  
  # Three components of the analytical log-posterior
  term1 <- log(sigma2) - log(denom) 
  term2 <- (sigma_mu_j_final * S_WR^2) / (sigma2 * denom)
  term3 <- - S_WR2 / sigma2  # R^2 term is mandatory due to inter-leaf C_i variance
  
  log_post <- 0.5 * sum(term1 + term2 + term3)
  
  return(log_post)
}





#' Direct Laplace Approximate Sampling for Leaf Node Parameters
#' 
#' @description Updates the 'mu' parameters of the terminal nodes by 
#' sampling from a Gaussian approximation centered at the posterior mode.
simulate_mu_skew_laplace <- function(tree, R, lambda1, gamma2, sigma2, sigma2_mu, 
                                     common_vars, aux_factor_var) {
  
  # 1. Identify terminal nodes
  which_terminal = as.numeric(which(tree$tree_matrix[,'terminal'] == 1))
  
  # 2. Identify terminals that violate ancestry constraints
  if(nrow(tree$tree_matrix) != 1) {
    terminal_ancestors = get_ancestors(tree)
    aux_table = table(terminal_ancestors[,1], terminal_ancestors[,2])
    which_terminal_no_double_split = NULL
    for (k in 1:nrow(aux_table)){
      split_var_ancestors = names(aux_table[k,])[aux_table[k,] != 0]
      invalid_ancestors = any(unlist(lapply(aux_factor_var, function(x) all(split_var_ancestors %in% x))))
      if (invalid_ancestors == TRUE) {
        which_terminal_no_double_split[k] = rownames(aux_table)[k]
      }
    }
    which_terminal_no_double_split = as.numeric(which_terminal_no_double_split[!is.na(which_terminal_no_double_split)]) 
  } else {
    which_terminal_no_double_split = 0
  }
  
  sigma2_mu_aux = rep(sigma2_mu, length(which_terminal))
  sigma2_mu_aux[which(which_terminal %in% which_terminal_no_double_split)] = 0
  
  node_sizes <- as.numeric(tree$tree_matrix[which_terminal, 'node_size'])
  mu_values  <- numeric(length(node_sizes))
  
  unique_leaf_indices <- sort(unique(tree$node_indices))
  gamma2_safe <- max(gamma2, 1e-10) # Numerical safety for skewness parameter
  
  # 3. Loop through each leaf to sample mu
  for(i in 1:length(node_sizes)) {
    
    if (sigma2_mu_aux[i] == 0) {
      mu_values[i] <- 0
      next
    }
    
    # Extract residuals and weights belonging to this leaf
    ind <- which(tree$node_indices == unique_leaf_indices[i])
    R_leaf <- R[ind]
    lambda_leaf <- lambda1[ind]
    
    # Numerical cleaning
    valid_mask <- is.finite(R_leaf) & is.finite(lambda_leaf)
    R_leaf <- R_leaf[valid_mask]
    lambda_leaf <- lambda_leaf[valid_mask]
    
    if (length(R_leaf) == 0) {
      mu_values[i] <- rnorm(1, 0, sqrt(sigma2_mu))
      next
    }
    
    # --- Mode Finding via IRLS (Iteratively Reweighted Least Squares) ---
    
    # Cold Start: Simple weighted average
    sum_lam <- sum(lambda_leaf)
    if(sum_lam <= 0) sum_lam <- 1e-10
    mu_hat <- sum(lambda_leaf * R_leaf) / (sum_lam + sigma2 / sigma2_mu)
    
    # Adaptive IRLS loop
    tolerance <- 1e-4
    max_iter_irls <- 10
    for(iter in 1:max_iter_irls) { 
      mu_old <- mu_hat
      C_i_tmp <- ifelse(R_leaf > mu_hat, 1/gamma2_safe, gamma2_safe)
      W_i <- C_i_tmp * lambda_leaf
      mu_hat <- sum(W_i * R_leaf) / (sum(W_i) + sigma2 / sigma2_mu)
      
      if(abs(mu_hat - mu_old) < tolerance) break
    }
    
    if (!is.finite(mu_hat)) mu_hat <- median(R_leaf)
    
    # --- Construct Laplace Approximation Parameters ---
    
    C_i_final <- ifelse(R_leaf > mu_hat, 1/gamma2_safe, gamma2_safe)
    W_i_final <- C_i_final * lambda_leaf
    S_W <- sum(W_i_final)
    
    laplace_mean <- mu_hat  
    laplace_var  <- sigma2 / (S_W + sigma2 / sigma2_mu)
    
    # Safety checks for variance
    if(!is.finite(laplace_var) || laplace_var <= 1e-14) laplace_var <- 1e-14
    laplace_sd <- sqrt(laplace_var)
    
    # --- Sample from Approximation ---
    # High-quality sample obtained via Gaussian approximation at the mode
    mu_star <- rnorm(1, mean = laplace_mean, sd = laplace_sd)
    
    # Final sanity check
    if(!is.finite(mu_star)) mu_star <- laplace_mean 
    
    mu_values[i] <- mu_star
  }
  
  # 4. Map simulated mu values back to the tree matrix
  tree$tree_matrix[,'mu'] = NA
  tree$tree_matrix[which_terminal, 'mu'] = mu_values
  tree$tree_matrix[which_terminal_no_double_split, 'mu'] = 0 
  
  return(tree)
}




