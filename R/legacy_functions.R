# Legacy Pure-R Implementation of SkT-BART
#
# This file contains a transparent pure-R implementation of SkT-BART
# retained for reproducibility and runtime comparison.
#
# The implementation uses BART-style local tree moves (grow, prune,
# and change) and a Laplace-based proposal correction for terminal-node
# parameters.
#
# For routine use, the Rcpp-accelerated implementation provided by
# `sktbart()` is recommended.


# ------------------------------------------------------------------------------
# Main legacy SkT-BART sampler
#
# Pure-R implementation used in the early development stage.
# Stores posterior samples of trees, error-distribution parameters,
# fitted values, and variable inclusion statistics.
#
# Recommended only for reproducibility and runtime comparisons.
# ------------------------------------------------------------------------------

sktbart_legacy <- function(
    y, x, sparse = FALSE, ntrees = 50, node_min_size = 5, gamma2 = 1,
    alpha = 0.95, beta = 2, lambda = NULL, d = 0.1, v = 30, mu_mu = 0,
    sigma2_mu = NULL, sigma2 = 1, b = 0.5, nburn = 1000, npost = 1000,
    nthin = 1, sigma_alpha = 1.5
) {
  if (is.null(sigma2_mu)) sigma2_mu <- 0.5^2 / ntrees
  if (is.null(lambda)) lambda <- rep(1, length(y))
  
  aux_identify_factor_variables <- NULL
  common_variables <- NULL
  
  total_iterations <- nburn + npost * nthin
  n <- length(y)
  p <- ncol(x)
  
  tree_store <- vector("list", npost)
  gamma2_store <- rep(NA_real_, npost)
  sigma2_store <- rep(NA_real_, npost)
  y_hat_store <- matrix(NA_real_, ncol = n, nrow = npost)
  bart_store <- matrix(NA_real_, ncol = n, nrow = npost)
  var_count <- rep(0, p)
  var_count_store <- matrix(0, ncol = p, nrow = npost)
  s_prob_store <- matrix(0, ncol = p, nrow = npost)
  v_store <- rep(NA_real_, npost)
  lambda_store <- matrix(NA_real_, ncol = n, nrow = npost)
  tree_accept_store <- rep(NA_real_, npost)
  
  tree_fits_cache <- matrix(0, ncol = ntrees, nrow = n)
  
  y_mean <- mean(y)
  y_sd <- sd(y)
  y_scale <- if (y_sd > 1e-8) (y - y_mean) / y_sd else (y - y_mean)
  
  lm_fit <- lm(y_scale ~ as.matrix(x))
  sigma_hat_sq <- var(lm_fit$residuals)
  sigma_beta <- 0.292 * sigma_hat_sq
  
  s <- rep(1 / p, p)
  curr_trees <- create_stump(num_trees = ntrees, y = y_scale, X = x)
  yhat_bart <- get_predictions(curr_trees, x, single_tree = (ntrees == 1))
  
  pb <- utils::txtProgressBar(min = 1, max = total_iterations, style = 3, width = 60)
  
  for (i in seq_len(total_iterations)) {
    utils::setTxtProgressBar(pb, i)
    accepted_this_iter <- 0L
    
    for (j in seq_len(ntrees)) {
      current_partial_residuals <- y_scale - yhat_bart + tree_fits_cache[, j]
      
      upd <- joint_update_tree_leaf_legacy(
        curr_tree = curr_trees[[j]],
        R = current_partial_residuals,
        X = x,
        lambda = lambda,
        gamma2 = gamma2,
        sigma2 = sigma2,
        sigma2_mu = sigma2_mu,
        alpha = alpha,
        beta = beta,
        node_min_size = node_min_size,
        s = s,
        iter = i,
        nburn = nburn,
        common_vars = common_variables,
        aux_factor_var = aux_identify_factor_variables
      )
      
      curr_trees[[j]] <- upd$tree
      accepted_this_iter <- accepted_this_iter + as.integer(upd$accepted)
      
      if (upd$accepted) {
        if (upd$move_type == "change" && length(curr_trees[[j]]$var) >= 2) {
          new_var <- curr_trees[[j]]$var[1]
          old_var <- curr_trees[[j]]$var[2]
          if (new_var %in% seq_len(p)) var_count[new_var] <- var_count[new_var] + 1
          if (old_var %in% seq_len(p)) var_count[old_var] <- var_count[old_var] - 1
        }
        if (upd$move_type == "grow" && length(curr_trees[[j]]$var) > 0) {
          for (vv in curr_trees[[j]]$var) {
            if (vv %in% seq_len(p)) var_count[vv] <- var_count[vv] + 1
          }
        }
        if (upd$move_type == "prune" && length(curr_trees[[j]]$var) > 0) {
          for (vv in curr_trees[[j]]$var) {
            if (vv %in% seq_len(p)) var_count[vv] <- var_count[vv] - 1
          }
        }
      }
      
      current_fit <- get_predictions(curr_trees[j], x, single_tree = TRUE)
      yhat_bart <- yhat_bart - tree_fits_cache[, j] + current_fit
      tree_fits_cache[, j] <- current_fit
    }
    
    current_residuals <- y_scale - yhat_bart
    
    error_params <- skewt_parameter_update(
      current_residuals,
      gamma2 = gamma2,
      v = v,
      d = d,
      lambda = lambda,
      sigma_alpha = sigma_alpha,
      sigma_beta = sigma_beta
    )
    
    v <- error_params$v
    sigma2 <- error_params$sigma2
    lambda <- error_params$lambda
    gamma2 <- error_params$gamma2
    
    if (isTRUE(sparse) && i > floor(total_iterations * 0.1)) {
      s <- update_s_(var_count, p, 1)
    }
    
    if ((i > nburn) && ((i - nburn) %% nthin == 0)) {
      idx <- (i - nburn) / nthin
      tree_store[[idx]] <- curr_trees
      sigma2_store[idx] <- sigma2
      y_hat_store[idx, ] <- yhat_bart
      bart_store[idx, ] <- yhat_bart
      var_count_store[idx, ] <- var_count
      s_prob_store[idx, ] <- s
      lambda_store[idx, ] <- lambda
      v_store[idx] <- v
      gamma2_store[idx] <- gamma2
      tree_accept_store[idx] <- accepted_this_iter / ntrees
    }
  }
  
  close(pb)
  cat("\n MCMC Sampling Complete. \n")
  
  results <- list(
    trees = tree_store,
    sigma2 = sigma2_store * y_sd^2,
    y_hat = y_hat_store * y_sd + y_mean,
    bart_hat = bart_store * y_sd + y_mean,
    v = v_store,
    gamma2 = gamma2_store,
    npost = npost,
    nburn = nburn,
    ntrees = ntrees,
    y_mean = y_mean,
    y_sd = y_sd,
    var_count_store = var_count_store,
    variable_weights = s_prob_store,
    lambda = lambda_store,
    tree_accept_rate = tree_accept_store
  )
  
  class(results) <- "sktbart"
  results
}


# ------------------------------------------------------------------------------
# Identify invalid terminal nodes under factor-splitting constraints
#
# This helper checks whether any terminal node violates the predefined
# factor-variable splitting rule. For each terminal node, it collects the
# splitting variables used along its ancestor path and compares them with
# `aux_factor_var`.
#
# Returns a numeric vector of terminal-node IDs that should be treated as
# invalid. If no constraint is supplied, or if the tree is a stump, it returns
# an empty numeric vector.
# ------------------------------------------------------------------------------
get_invalid_leaf_ids_legacy <- function(tree, common_vars = NULL, aux_factor_var = NULL) {
  terminal_nodes <- as.numeric(which(tree$tree_matrix[, "terminal"] == 1))
  if (nrow(tree$tree_matrix) == 1 || is.null(aux_factor_var) || length(aux_factor_var) == 0) {
    return(numeric(0))
  }
  terminal_ancestors <- get_ancestors(tree)
  if (is.null(terminal_ancestors) || nrow(terminal_ancestors) == 0) {
    return(numeric(0))
  }
  
  ancestor_table <- table(terminal_ancestors[, 1], terminal_ancestors[, 2])
  invalid <- numeric(0)
  for (k in seq_len(nrow(ancestor_table))) {
    ancestor_vars <- names(ancestor_table[k, ])[ancestor_table[k, ] != 0]
    is_invalid <- any(vapply(aux_factor_var, function(z) all(ancestor_vars %in% z), logical(1)))
    if (is_invalid) invalid <- c(invalid, as.numeric(rownames(ancestor_table)[k]))
  }
  intersect(terminal_nodes, invalid)
}

                             
# Joint tree-leaf update (legacy pure-R version)
#
# A candidate tree structure is generated using standard BART local moves.
# Conditional on the proposed tree, terminal-node parameters are proposed
# from a Laplace-based Gaussian proposal distribution.
#
# The acceptance ratio evaluates the joint skew-t posterior kernel and
# includes the forward/reverse Laplace proposal-density correction for
# terminal-node parameters.
                             
joint_update_tree_leaf_legacy <- function(curr_tree, R, X, lambda, gamma2, sigma2,
                                          sigma2_mu, alpha, beta, node_min_size, s,
                                          iter, nburn,
                                          common_vars = NULL, aux_factor_var = NULL) {
  move_type <- sample_move(curr_tree, iter, nburn)
  
  proposed_structure <- update_tree(
    y = R,
    X = X,
    type = move_type,
    curr_tree = curr_tree,
    node_min_size = node_min_size,
    s = s,
    common_vars = common_vars,
    aux_factor_var = aux_factor_var
  )
  
  old_prop <- leaf_laplace_proposal_legacy(
    tree = curr_tree,
    R = R,
    lambda = lambda,
    gamma2 = gamma2,
    sigma2 = sigma2,
    sigma2_mu = sigma2_mu,
    common_vars = common_vars,
    aux_factor_var = aux_factor_var,
    draw = FALSE,
    eval_tree = curr_tree
  )
  
  new_prop <- leaf_laplace_proposal_legacy(
    tree = proposed_structure,
    R = R,
    lambda = lambda,
    gamma2 = gamma2,
    sigma2 = sigma2,
    sigma2_mu = sigma2_mu,
    common_vars = common_vars,
    aux_factor_var = aux_factor_var,
    draw = TRUE
  )
  proposed_tree <- new_prop$tree
  
  log_target_old <- log_joint_tree_leaf_skewt_legacy(
    tree = curr_tree,
    R = R,
    lambda = lambda,
    gamma2 = gamma2,
    sigma2 = sigma2,
    sigma2_mu = sigma2_mu,
    alpha = alpha,
    beta = beta,
    common_vars = common_vars,
    aux_factor_var = aux_factor_var
  )
  
  log_target_new <- log_joint_tree_leaf_skewt_legacy(
    tree = proposed_tree,
    R = R,
    lambda = lambda,
    gamma2 = gamma2,
    sigma2 = sigma2,
    sigma2_mu = sigma2_mu,
    alpha = alpha,
    beta = beta,
    common_vars = common_vars,
    aux_factor_var = aux_factor_var
  )
  log_ratio <- log_target_new - log_target_old +
    old_prop$log_q - new_prop$log_q
  
  accepted <- FALSE
  out_tree <- curr_tree
  if (!isTRUE(proposed_tree$ForceStump) &&
      is.finite(log_ratio) &&
      (log_ratio >= 0 || log(runif(1)) < log_ratio)) {
    out_tree <- proposed_tree
    accepted <- TRUE
  }
  
  list(tree = out_tree, accepted = accepted, move_type = move_type, log_ratio = log_ratio)
}

                            
# ------------------------------------------------------------------------------
# Laplace proposal construction for terminal-node parameters
#
# For each terminal node, a local mode is obtained using an
# iteratively reweighted least squares (IRLS) approximation.
# A Gaussian proposal is then constructed using the local curvature.
#
# This proposal is used only for generating candidate leaf parameters.
# ------------------------------------------------------------------------------
                             
leaf_laplace_proposal_legacy <- function(tree, R, lambda, gamma2, sigma2, sigma2_mu,
                                         common_vars = NULL, aux_factor_var = NULL,
                                         draw = TRUE, eval_tree = NULL) {
  if (is.null(eval_tree)) eval_tree <- tree
  
  leaf_ids <- sort(unique(tree$node_indices))
  invalid_leaf_ids <- get_invalid_leaf_ids_legacy(tree, common_vars, aux_factor_var)
  gamma2_safe <- max(gamma2, 1e-10)
  
  out_tree <- tree
  out_tree$tree_matrix[, "mu"] <- NA
  log_q <- 0
  
  for (leaf_id in leaf_ids) {
    obs <- which(tree$node_indices == leaf_id)
    
    if (leaf_id %in% invalid_leaf_ids) {
      mu_value <- 0
      out_tree$tree_matrix[leaf_id, "mu"] <- 0
      next
    }
    
    if (length(obs) == 0) {
      mu_hat <- 0
      laplace_var <- sigma2_mu
    } else {
      R_leaf <- R[obs]
      lambda_leaf <- lambda[obs]
      ok <- is.finite(R_leaf) & is.finite(lambda_leaf)
      R_leaf <- R_leaf[ok]
      lambda_leaf <- lambda_leaf[ok]
      
      if (length(R_leaf) == 0) {
        mu_hat <- 0
        laplace_var <- sigma2_mu
      } else {
        mu_hat <- sum(lambda_leaf * R_leaf) /
          (sum(lambda_leaf) + sigma2 / sigma2_mu)
        
        for (iter in seq_len(10)) {
          mu_old <- mu_hat
          w <- ifelse(R_leaf > mu_hat, 1 / gamma2_safe, gamma2_safe) * lambda_leaf
          mu_hat <- sum(w * R_leaf) / (sum(w) + sigma2 / sigma2_mu)
          if (!is.finite(mu_hat)) mu_hat <- mu_old
          if (abs(mu_hat - mu_old) < 1e-4) break
        }
        
        w_final <- ifelse(R_leaf > mu_hat, 1 / gamma2_safe, gamma2_safe) * lambda_leaf
        laplace_var <- sigma2 / (sum(w_final) + sigma2 / sigma2_mu)
      }
    }
    
    laplace_sd <- sqrt(max(laplace_var, 1e-14))
    if (draw) {
      mu_value <- rnorm(1, mean = mu_hat, sd = laplace_sd)
      if (!is.finite(mu_value)) mu_value <- mu_hat
    } else {
      mu_value <- eval_tree$tree_matrix[leaf_id, "mu"]
      if (!is.finite(mu_value)) mu_value <- mu_hat
    }
    
    out_tree$tree_matrix[leaf_id, "mu"] <- mu_value
    log_q <- log_q + dnorm(mu_value, mean = mu_hat, sd = laplace_sd, log = TRUE)
  }
  
  list(tree = out_tree, log_q = log_q)
}


# ------------------------------------------------------------------------------
# Joint posterior kernel evaluation
#
# Computes the unnormalized joint posterior kernel for a tree and its
# terminal-node parameters under the skew-t likelihood and Gaussian
# leaf prior.
#
# Used in the Metropolis-Hastings acceptance ratio.
# ------------------------------------------------------------------------------

log_joint_tree_leaf_skewt_legacy <- function(tree, R, lambda, gamma2, sigma2, sigma2_mu,
                                             alpha, beta,
                                             common_vars = NULL, aux_factor_var = NULL) {
  terminal_nodes <- as.numeric(which(tree$tree_matrix[, "terminal"] == 1))
  invalid_leaf_ids <- get_invalid_leaf_ids_legacy(tree, common_vars, aux_factor_var)
  gamma2_safe <- max(gamma2, 1e-10)
  
  log_lik <- 0
  log_mu_prior <- 0
  
  for (leaf_id in terminal_nodes) {
    mu_leaf <- tree$tree_matrix[leaf_id, "mu"]
    if (!is.finite(mu_leaf)) return(-Inf)
    
    if (leaf_id %in% invalid_leaf_ids) {
      if (abs(mu_leaf) > 1e-12) return(-Inf)
    } else {
      log_mu_prior <- log_mu_prior + dnorm(mu_leaf, 0, sqrt(sigma2_mu), log = TRUE)
    }
    
    obs <- which(tree$node_indices == leaf_id)
    if (length(obs) > 0) {
      resid <- R[obs] - mu_leaf
      c_i <- ifelse(R[obs] >= mu_leaf, 1 / gamma2_safe, gamma2_safe)
      log_lik <- log_lik - 0.5 * sum(lambda[obs] * c_i * resid^2) / sigma2
    }
  }
  
  log_lik + log_mu_prior + get_tree_prior(tree, alpha, beta, common_vars)
}

                             
# ------------------------------------------------------------------------------
# Legacy leaf-parameter updater
#
# Earlier Laplace-based leaf update routine retained for reference.
# Superseded by the joint tree-leaf proposal mechanism implemented in
# `joint_update_tree_leaf_legacy()`.
# ------------------------------------------------------------------------------

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

fill_tree_details = function(curr_tree, X) {
  
  # Collect right bits of tree
  tree_matrix = curr_tree$tree_matrix
  
  # Create a new tree matrix to overwrite
  new_tree_matrix = tree_matrix
  
  # Start with dummy node indices
  node_indices = rep(1, nrow(X))
  
  # For all but the top row, find the number of observations falling into each one
  for(i in 2:nrow(tree_matrix)) {
    
    # Get the parent
    curr_parent = as.numeric(tree_matrix[i,'parent'])
    
    # Find the split variable and value of the parent
    split_var = as.numeric(tree_matrix[curr_parent,'split_variable'])
    split_val = as.numeric(tree_matrix[curr_parent, 'split_value'])
    
    # Find whether it's a left or right terminal node
    left_or_right = ifelse(tree_matrix[curr_parent,'child_left'] == i,
                           'left', 'right')
    if(left_or_right == 'left') {
      # If left use less than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] < split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] < split_val] = i
    } else {
      # If right use greater than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] >= split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] >= split_val] = i
    }
  } # End of loop through table
  
  return(list(tree_matrix = new_tree_matrix,
              node_indices = node_indices))
  
} # End of function
