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


update_v_skewt <- function(n, d, lambda1){
  eta = sum(lambda1 - log(lambda1))/2+d
  T=n
  f <- function(x) {
    (T / 2) * (log(x / 2) + 1 - digamma(x / 2)) + (1 / x) - eta
  }
  f_lower = f(1e-10)
  f_upper = f(100)
  
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

tree_full_skewt_laplace <- function(tree, R, lambda, sigma2, sigma2_mu, gamma2, 
                                    common_vars, aux_factor_var) {
  
  # === 留在 R 中的部分 (处理繁杂的拓扑验证) ===
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
  unique_leaf_ids <- terminal_nodes
  sigma2_mu_j <- rep(sigma2_mu, length(unique_leaf_ids))
  
  if (length(terminals_invalid) > 0 && any(terminals_invalid != 0)) {
    invalid_idx <- which(unique_leaf_ids %in% terminals_invalid)
    sigma2_mu_j[invalid_idx] <- 0
  }
  
  leaf_ids <- as.integer(tree$node_indices)
  
  # === 抛给 C++ 的部分 (处理高密度数学运算) ===
  # Step 4-6: C++ backend
  marg_lik <- laplace_irls_cpp(R = R,
                               lambda = lambda,
                               leaf_ids = leaf_ids,
                               unique_leaf_ids = as.integer(unique_leaf_ids),
                               sigma2_mu_j = sigma2_mu_j,
                               gamma2 = gamma2,
                               sigma2 = sigma2)
  
  return(marg_lik)
}







simulate_mu_skew_laplace <- function(tree, R, lambda1, gamma2, sigma2, sigma2_mu, 
                                     common_vars, aux_factor_var) {
  
  # 1. 在 R 层面处理繁杂的树节点拓扑和多重分裂验证逻辑 (保持原样)
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
  
  unique_leaf_indices <- sort(unique(tree$node_indices))
  
  # 对齐参数尺寸：给所有现存的叶子节点分配先验方差
  sigma2_mu_aux <- rep(sigma2_mu, length(unique_leaf_indices))
  
  # 将违规节点的方差强行置为 0
  if (length(which_terminal_no_double_split) > 0 && any(which_terminal_no_double_split != 0)) {
    invalid_idx <- which(unique_leaf_indices %in% which_terminal_no_double_split)
    sigma2_mu_aux[invalid_idx] <- 0
  }
  
  # 2. 调用 C++ 后端：极限压缩 IRLS 抽样循环的时间
  mu_values <- simulate_mu_irls_cpp(
    R = R,
    lambda = lambda1,
    leaf_ids = as.integer(tree$node_indices),
    unique_leaf_ids = as.integer(unique_leaf_indices),
    sigma2_mu_j = sigma2_mu_aux,
    gamma2 = gamma2,
    sigma2 = sigma2
  )
  
  # 3. 把算好的参数原封不动地写回树矩阵里
  tree$tree_matrix[, 'mu'] <- NA
  tree$tree_matrix[unique_leaf_indices, 'mu'] <- mu_values
  
  # 防御性覆写，确保无效节点强行归零
  if (any(which_terminal_no_double_split != 0)) {
    tree$tree_matrix[which_terminal_no_double_split, 'mu'] <- 0 
  }
  
  return(tree)
}
