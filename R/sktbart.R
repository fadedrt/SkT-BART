
sktbart <- function(
    y,
    x,
    sparse        = FALSE,
    ntrees        = 50,
    node_min_size = 5,
    gamma2        = 1,
    alpha         = 0.95,
    beta          = 2,
    lambda        = NULL,
    d             = 0.1,
    v             = 30,
    mu_mu         = 0,
    sigma2_mu     = NULL,
    sigma2        = 1,
    b             = 0.5,
    nburn         = 1000,
    npost         = 1000,
    nthin         = 1,
    sigma_alpha   = 1.5
) {
  
  # --- 1. Hyperparameter Initialization ---
  if (is.null(sigma2_mu)) {
    sigma2_mu <- 0.5^2 / ntrees  
  }
  if (is.null(lambda)) {
    lambda <- rep(1, length(y))
  }
  
  # Placeholders for constraints and factor handling
  aux_identify_factor_variables <- NULL
  common_variables              <- NULL
  
  total_iterations <- nburn + npost * nthin
  n <- length(y)
  p <- ncol(x)
  
  # --- 2. Storage Allocation ---
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
  
  # Local cache for tree fits to speed up residual calculations
  tree_fits_cache <- matrix(0, ncol = ntrees, nrow = n)
  
  # --- 3. Data Pre-processing ---
  y_mean <- mean(y)
  y_sd   <- sd(y) 
  
  # Scale response: handle constant y case
  y_scale <- if(y_sd > 1e-8) (y - y_mean) / y_sd else (y - y_mean)
  
  # Adaptive prior for sigma2 based on linear model fit
  lm_fit       <- lm(y_scale ~ as.matrix(x))
  sigma_hat_sq <- var(lm_fit$residuals)
  sigma_beta   <- 0.292 * sigma_hat_sq  
  
  # Initialize variable selection probabilities
  s <- rep(1 / p, p)
  
  # --- 4. Model Initialization ---
  curr_trees <- create_stump(num_trees = ntrees, y = y_scale, X = x)
  new_trees  <- curr_trees
  yhat_bart  <- get_predictions(curr_trees, x, single_tree = (ntrees == 1))
  
  # Initialize progress bar
  pb <- utils::txtProgressBar(min = 1, max = total_iterations, style = 3, width = 60)
  
  # --- 5. Main MCMC Loop ---
  for (i in seq_len(total_iterations)) {
    utils::setTxtProgressBar(pb, i)
    
    # 5.1 Tree Ensemble Update (Back-fitting)
    for (j in seq_len(ntrees)) {
      
      # Partial residuals (excluding current tree)
      current_partial_residuals <- y_scale - yhat_bart + tree_fits_cache[, j]
      
      # Propose structure move
      move_type <- sample_move(curr_trees[[j]], i, nburn)
      new_trees[[j]] <- update_tree(
        y              = current_partial_residuals, 
        X              = x,
        type           = move_type, 
        curr_tree      = curr_trees[[j]],
        node_min_size  = node_min_size, 
        s              = s,
        common_vars    = common_variables,
        aux_factor_var = aux_identify_factor_variables
      )
      
      # Laplace Marginal Likelihood Evaluation
      l_old <- tree_full_skewt_laplace(
        tree           = curr_trees[[j]], 
        R              = current_partial_residuals,
        lambda         = lambda, 
        sigma2         = sigma2, 
        sigma2_mu      = sigma2_mu, 
        gamma2         = gamma2,
        common_vars    = common_variables, 
        aux_factor_var = aux_identify_factor_variables
      ) + get_tree_prior(curr_trees[[j]], alpha, beta, common_variables)
      
      l_new <- tree_full_skewt_laplace(
        tree           = new_trees[[j]], 
        R              = current_partial_residuals,
        lambda         = lambda, 
        sigma2         = sigma2, 
        sigma2_mu      = sigma2_mu, 
        gamma2         = gamma2,
        common_vars    = common_variables, 
        aux_factor_var = aux_identify_factor_variables
      ) + get_tree_prior(new_trees[[j]], alpha, beta, common_variables)
      
      # Metropolis-Hastings Acceptance
      ratio <- if (isTRUE(new_trees[[j]]$ForceStump)) -Inf else l_new - l_old
      
      if (ratio > 0 || ratio > log(runif(1))) {
        curr_trees[[j]] <- new_trees[[j]]
        # Track variable importance
        if (move_type == 'change') {
          var_count[curr_trees[[j]]$var[1]] <- var_count[curr_trees[[j]]$var[1]] + 1
          var_count[curr_trees[[j]]$var[2]] <- var_count[curr_trees[[j]]$var[2]] - 1
        }
        if (move_type == 'grow')  var_count[curr_trees[[j]]$var] <- var_count[curr_trees[[j]]$var] + 1
        if (move_type == 'prune') var_count[curr_trees[[j]]$var] <- var_count[curr_trees[[j]]$var] - 1
      }
      
      # Update Leaf Node Parameters (Laplace Approximation)
      curr_trees[[j]] <- simulate_mu_skew_laplace(
        tree           = curr_trees[[j]], 
        R              = current_partial_residuals,
        lambda1        = lambda, 
        gamma2         = gamma2, 
        sigma2         = sigma2, 
        sigma2_mu      = sigma2_mu,
        common_vars    = common_variables, 
        aux_factor_var = aux_identify_factor_variables
      )
      
      # Sync total predictions
      current_fit <- get_predictions(curr_trees[j], x, single_tree = TRUE)
      yhat_bart <- yhat_bart - tree_fits_cache[, j] + current_fit
      tree_fits_cache[, j] <- current_fit
    }
    
    # 5.2 Skew-t Distribution Parameter Updates
    current_residuals <- y_scale - yhat_bart
    
    error_params <- skewt_parameter_update(
      current_residuals, 
      gamma2      = gamma2, 
      v           = v, 
      d           = d, 
      lambda      = lambda, 
      sigma_alpha = sigma_alpha, 
      sigma_beta  = sigma_beta
    )
    
    v      <- error_params$v
    sigma2 <- error_params$sigma2
    lambda <- error_params$lambda
    gamma2 <- error_params$gamma2
    
    # 5.3 Dirichlet Sparsity Update
    if(isTRUE(sparse) && i > floor(total_iterations * 0.1)){
      s <- update_s(var_count, p, 1)
    }
    
    # 5.4 Posterior Collection
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
  
  # --- 6. Results Formatting (Rescale) ---
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
