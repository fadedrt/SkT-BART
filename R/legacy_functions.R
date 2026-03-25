#' @title Skew-t Bayesian Additive Regression Trees (SkTBART)
#' @description 
#' Complete Pure R implementation of the legacy SkTBART algorithm. 
#' NOTE: This version does NOT use C++ acceleration (Pure R implementation).
#' It handles non-conjugate Skew-t error distributions via Laplace approximations 
#' to robustly model data with outliers and asymmetric residuals.
# ==============================================================================
# 1. LAPLACE APPROXIMATION LIKELIHOOD & SAMPLING
# ==============================================================================

#' @export
tree_full_skewt_laplace_legacy <- function(tree, R, lambda, sigma2, sigma2_mu, gamma2, 
                                           common_vars, aux_factor_var) {
  
  terminal_nodes <- as.numeric(which(tree$tree_matrix[, 'terminal'] == 1))
  
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
  
  sigma2_mu_named <- rep(sigma2_mu, length(terminal_nodes))
  names(sigma2_mu_named) <- as.character(terminal_nodes)
  if (length(terminals_invalid) > 0 && any(terminals_invalid != 0)) {
    sigma2_mu_named[as.character(terminals_invalid)] <- 0
  }
  
  leaf_ids <- as.character(tree$node_indices)
  
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
  
  C_i <- ifelse(R > mu_hat, 1 / gamma2, gamma2)
  W_i <- C_i * lambda
  
  S_W   <- tapply(W_i, leaf_ids, sum)
  S_WR  <- tapply(W_i * R, leaf_ids, sum)
  S_WR2 <- tapply(W_i * R^2, leaf_ids, sum)
  
  sigma_mu_j_final <- sigma2_mu_named[names(S_W)]
  
  denom <- S_W * sigma_mu_j_final + sigma2
  term1 <- log(sigma2) - log(denom) 
  term2 <- (sigma_mu_j_final * S_WR^2) / (sigma2 * denom)
  term3 <- - S_WR2 / sigma2  
  
  return(0.5 * sum(term1 + term2 + term3))
}

#' @export
simulate_mu_skew_laplace_legacy <- function(tree, R, lambda1, gamma2, sigma2, sigma2_mu, 
                                            common_vars, aux_factor_var) {
  
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
    
    W_final <- ifelse(R_leaf > mu_hat, 1/gamma2_safe, gamma2_safe) * lambda_leaf
    laplace_var <- sigma2 / (sum(W_final) + sigma2 / sigma2_mu)
    laplace_sd <- sqrt(max(laplace_var, 1e-14))
    
    mu_star <- rnorm(1, mean = mu_hat, sd = laplace_sd)
    mu_values[i] <- ifelse(is.finite(mu_star), mu_star, mu_hat)
  }
  
  tree$tree_matrix[, 'mu'] <- NA
  tree$tree_matrix[which_terminal, 'mu'] <- mu_values
  tree$tree_matrix[which_terminal_no_double_split, 'mu'] <- 0 
  
  return(tree)
}

# ==============================================================================
# 2. MAIN MODEL FITTING FUNCTION
# ==============================================================================

#' @export
sktbart_legacy <- function(
    y, x, sparse = FALSE, ntrees = 50, node_min_size = 5, gamma2 = 1,
    alpha = 0.95, beta = 2, lambda = NULL, d = 0.1, v = 30, mu_mu = 0,
    sigma2_mu = NULL, sigma2 = 1, b = 0.5, nburn = 1000, npost = 1000,
    nthin = 1, sigma_alpha = 1.5
) {
  
  if (is.null(sigma2_mu)) sigma2_mu <- 0.5^2 / ntrees  
  if (is.null(lambda)) lambda <- rep(1, length(y))
  
  aux_identify_factor_variables <- NULL
  common_variables              <- NULL
  
  total_iterations <- nburn + npost * nthin
  n <- length(y)
  p <- ncol(x)
  
  tree_store      <- vector('list', npost)
  gamma2_store    <- rep(NA, npost)
  sigma2_store    <- rep(NA, npost)
  y_hat_store     <- matrix(NA, ncol = n, nrow = npost)
  bart_store      <- matrix(NA, ncol = n, nrow = npost)
  var_count       <- rep(0, p)
  var_count_store <- matrix(0, ncol = p, nrow = npost)
  s_prob_store    <- matrix(0, ncol = p, nrow = npost)
  v_store         <- rep(NA, npost)
  lambda_store    <- matrix(NA, ncol = n, nrow = npost)
  
  tree_fits_cache <- matrix(0, ncol = ntrees, nrow = n)
  
  y_mean <- mean(y)
  y_sd   <- sd(y) 
  y_scale <- if(y_sd > 1e-8) (y - y_mean) / y_sd else (y - y_mean)
  
  lm_fit       <- lm(y_scale ~ as.matrix(x))
  sigma_hat_sq <- var(lm_fit$residuals)
  sigma_beta   <- 0.292 * sigma_hat_sq  
  
  s <- rep(1 / p, p)
  
  # Initialize stump (Assumes create_stump is available in your env)
  curr_trees <- create_stump(num_trees = ntrees, y = y_scale, X = x)
  new_trees  <- curr_trees
  yhat_bart  <- get_predictions(curr_trees, x, single_tree = (ntrees == 1))
  
  pb <- utils::txtProgressBar(min = 1, max = total_iterations, style = 3, width = 60)
  
  for (i in seq_len(total_iterations)) {
    utils::setTxtProgressBar(pb, i)
    
    for (j in seq_len(ntrees)) {
      current_partial_residuals <- y_scale - yhat_bart + tree_fits_cache[, j]
      
      move_type <- sample_move(curr_trees[[j]], i, nburn)
      
      # Assumes update_tree is available in your env
      new_trees[[j]] <- update_tree(
        y = current_partial_residuals, X = x, type = move_type, 
        curr_tree = curr_trees[[j]], node_min_size = node_min_size, 
        s = s, common_vars = common_variables, aux_factor_var = aux_identify_factor_variables
      )
      
      l_old <- tree_full_skewt_laplace_legacy(
        tree = curr_trees[[j]], R = current_partial_residuals, lambda = lambda, 
        sigma2 = sigma2, sigma2_mu = sigma2_mu, gamma2 = gamma2,
        common_vars = common_variables, aux_factor_var = aux_identify_factor_variables
      ) + get_tree_prior(curr_trees[[j]], alpha, beta, common_variables)
      
      l_new <- tree_full_skewt_laplace_legacy(
        tree = new_trees[[j]], R = current_partial_residuals, lambda = lambda, 
        sigma2 = sigma2, sigma2_mu = sigma2_mu, gamma2 = gamma2,
        common_vars = common_variables, aux_factor_var = aux_identify_factor_variables
      ) + get_tree_prior(new_trees[[j]], alpha, beta, common_variables)
      
      ratio <- if (isTRUE(new_trees[[j]]$ForceStump)) -Inf else l_new - l_old
      
      if (ratio > 0 || ratio > log(runif(1))) {
        curr_trees[[j]] <- new_trees[[j]]
        if (move_type == 'change') {
          var_count[curr_trees[[j]]$var[1]] <- var_count[curr_trees[[j]]$var[1]] + 1
          var_count[curr_trees[[j]]$var[2]] <- var_count[curr_trees[[j]]$var[2]] - 1
        }
        if (move_type == 'grow')  var_count[curr_trees[[j]]$var] <- var_count[curr_trees[[j]]$var] + 1
        if (move_type == 'prune') var_count[curr_trees[[j]]$var] <- var_count[curr_trees[[j]]$var] - 1
      }
      
      curr_trees[[j]] <- simulate_mu_skew_laplace_legacy(
        tree = curr_trees[[j]], R = current_partial_residuals, lambda1 = lambda, 
        gamma2 = gamma2, sigma2 = sigma2, sigma2_mu = sigma2_mu,
        common_vars = common_variables, aux_factor_var = aux_identify_factor_variables
      )
      
      current_fit <- get_predictions(curr_trees[j], x, single_tree = TRUE)
      yhat_bart <- yhat_bart - tree_fits_cache[, j] + current_fit
      tree_fits_cache[, j] <- current_fit
    }
    
    current_residuals <- y_scale - yhat_bart
    
    # Assumes skewt_parameter_update is available in your env
    error_params <- skewt_parameter_update(
      current_residuals, gamma2 = gamma2, v = v, d = d, 
      lambda = lambda, sigma_alpha = sigma_alpha, sigma_beta = sigma_beta
    )
    
    v      <- error_params$v
    sigma2 <- error_params$sigma2
    lambda <- error_params$lambda
    gamma2 <- error_params$gamma2
    
    if(isTRUE(sparse) && i > floor(total_iterations * 0.1)){
      s <- update_s_(var_count, p, 1)
    }
    
    if ((i > nburn) && ((i - nburn) %% nthin == 0)) {
      idx <- (i - nburn) / nthin
      tree_store[[idx]]      <- curr_trees
      sigma2_store[idx]      <- sigma2
      y_hat_store[idx, ]     <- yhat_bart
      bart_store[idx, ]      <- yhat_bart
      var_count_store[idx, ] <- var_count
      s_prob_store[idx, ]    <- s
      lambda_store[idx, ]    <- lambda
      v_store[idx]           <- v
      gamma2_store[idx]      <- gamma2
    }
  }
  
  close(pb)
  cat('\n MCMC Sampling Complete. \n')
  
  results <- list(
    trees            = tree_store,
    sigma2           = sigma2_store * y_sd^2,
    y_hat            = y_hat_store * y_sd + y_mean,
    bart_hat         = bart_store * y_sd + y_mean,
    v                = v_store,
    gamma2           = gamma2_store,
    npost            = npost,
    nburn            = nburn,
    ntrees           = ntrees,
    y_mean           = y_mean,
    y_sd             = y_sd,
    var_count_store  = var_count_store,
    variable_weights = s_prob_store,
    lambda           = lambda_store
  )
  
  class(results) <- "sktbart"
  return(results)
}
