#' BART with Skew-t Errors via Laplace Approximation
#'
#' @param y A numeric vector of the response variable.
#' @param x A data frame or matrix of predictors.
#' @param sparse Logical; if TRUE, uses a Dirichlet prior for sparse variable selection.
#' @param ntrees Number of trees in the ensemble (default is 50).
#' @param node_min_size Minimum number of observations required in a terminal node.
#' @param alpha Base hyperparameter for tree prior (controls depth).
#' @param beta Power hyperparameter for tree prior (controls depth).
#' @param nu Degrees of freedom for the Skew-t distribution (initial value).
#' @param lambda Latent weights for the scale-mixture representation (initial vector).
#' @param mu_mu Mean of the prior for leaf node parameters.
#' @param sigma2 Initial error variance.
#' @param sigma2_mu Prior variance for leaf node parameters.
#' @param nburn Number of burn-in iterations for MCMC.
#' @param npost Number of posterior samples to collect.
#' @param nthin Thinning interval for MCMC.
#'
#' @return An object of class 'sktbart' containing posterior samples of trees, 
#' sigma2, gamma2, v, and predicted values.
#' 
#' @importFrom stats rgamma rexp dnorm sd rchisq rnorm pnorm as.formula terms xtabs lm var median
#' @importFrom truncnorm rtruncnorm
#' @importFrom dbarts makeModelMatrixFromDataFrame
#' @importFrom revdbayes ru
#' @export
#'
#' @examples
#' # Example usage with simulated data
#' # n <- 100; p <- 5
#' # x <- matrix(runif(n * p), n, p)
#' # y <- sin(pi * x[,1] * x[,2]) + rnorm(n)
#' # model <- bart_skewt(y, x, nburn = 100, npost = 100)
#'
bart_skewt <- function(y,
                               x,
                               sparse = FALSE,
                               ntrees = 50,
                               node_min_size = 5,
                               gamma2 = 1,
                               alpha = 0.95,
                               beta = 2,
                               lambda = NULL,
                               d = 0.1,
                               v = 30,
                               mu_mu = 0,
                               sigma2_mu = NULL,
                               sigma2 = 1,
                               b = 0.5,
                               nburn = 1000,
                               npost = 1000,
                               nthin = 1,
                               sigma_alpha = 1.5
) {
  # --- Hyperparameter Initialization ---
  if (is.null(sigma2_mu)) {
    sigma2_mu = 0.5^2 / ntrees  # Scaled prior for leaf nodes
  }
  if (is.null(lambda)) {
    lambda = rep(1, length(y))
  }
  
  # Placeholders for factor variable handling
  aux_identify_factor_variables = NULL
  common_variables = NULL
  
  total_iterations = nburn + npost * nthin
  n = length(y)
  p = ncol(x)
  
  # --- Initialize Storage ---
  tree_store      = vector('list', npost)
  gamma2_store    = rep(NA, npost)
  sigma2_store    = rep(NA, npost)
  y_hat_store     = matrix(NA, ncol = n, nrow = npost)
  bart_store      = matrix(NA, ncol = n, nrow = npost)
  var_count       = rep(0, p)
  var_count_store = matrix(0, ncol = p, nrow = npost)
  s_prob_store    = matrix(0, ncol = p, nrow = npost)
  v_store         = rep(NA, npost)
  lambda_store    = matrix(NA, ncol = n, nrow = npost)
  
  # Local cache for tree fits to speed up residual calculations
  tree_fits_store = matrix(0, ncol = ntrees, nrow = n)
  
  # --- Standardize Response Variable ---
  y_mean  = mean(y)
  y_sd    = sd(y) 
  # Handle zero-variance edge case
  y_scale = if(y_sd > 0) (y - y_mean) / y_sd else (y - y_mean)
  
  # Initial estimate of error variance from a linear model
  lm_fit       = lm(y_scale ~ x)
  sigma_hat_sq = var(lm_fit$residuals)
  sigma_beta   = 0.292 * sigma_hat_sq  # Prior scale for sigma2
  
  # Initialize variable inclusion probabilities (Uniform)
  s = rep(1 / p, p)
  
  # --- Initialize Tree Structures ---
  curr_trees = create_stump(num_trees = ntrees, y = y_scale, X = x)
  new_trees  = curr_trees
  yhat_bart  = get_predictions(curr_trees, x, single_tree = (ntrees == 1))
  
  # --- Progress Bar ---
  pb = utils::txtProgressBar(min = 1, max = total_iterations, style = 3, width = 60)
  
  # --- Main MCMC Loop ---
  for (i in seq_len(total_iterations)) {
    utils::setTxtProgressBar(pb, i)
    
    # 1. Update Trees (Back-fitting)
    for (j in seq_len(ntrees)) {
      # Compute partial residuals for the j-th tree
      current_partial_residuals = y_scale - yhat_bart + tree_fits_store[, j]
      
      # Propose a tree move (Grow, Prune, Change)
      type = sample_move(curr_trees[[j]], i, nburn)
      new_trees[[j]] = update_tree(
        y = current_partial_residuals, X = x,
        type = type, curr_tree = curr_trees[[j]],
        node_min_size = node_min_size, s = s,
        common_vars = common_variables,
        aux_factor_var = aux_identify_factor_variables
      )
      
      # Compute log-likelihood with Skew-t and Tree Prior
      l_old = tree_full_skewt_laplace(
        tree = curr_trees[[j]], R = current_partial_residuals,
        lambda = lambda, sigma2 = sigma2, sigma2_mu = sigma2_mu, gamma2 = gamma2,
        common_variables, aux_identify_factor_variables
      ) + get_tree_prior(curr_trees[[j]], alpha = alpha, beta = beta, common_variables)
      
      l_new = tree_full_skewt_laplace(
        tree = new_trees[[j]], R = current_partial_residuals,
        lambda = lambda, sigma2 = sigma2, sigma2_mu = sigma2_mu, gamma2 = gamma2,
        common_variables, aux_identify_factor_variables
      ) + get_tree_prior(new_trees[[j]], alpha = alpha, beta = beta, common_variables)
      
      # Metropolis-Hastings Step
      acceptance_prob = if (isTRUE(new_trees[[j]]$ForceStump)) 0 else l_new - l_old
      if (acceptance_prob > 0 || acceptance_prob > log(runif(1))) {
        curr_trees[[j]] = new_trees[[j]]
        # Update variable importance counts
        if (type == 'change') {
          var_count[curr_trees[[j]]$var[1]] = var_count[curr_trees[[j]]$var[1]] + 1
          var_count[curr_trees[[j]]$var[2]] = var_count[curr_trees[[j]]$var[2]] - 1
        }
        if (type == 'grow')  var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] + 1
        if (type == 'prune') var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] - 1
      }
      
      # Update Leaf Node Parameters (mu)
      curr_trees[[j]] = simulate_mu_skew_laplace(
        tree = curr_trees[[j]], R = current_partial_residuals,
        lambda = lambda, gamma2 = gamma2, sigma2 = sigma2, sigma2_mu = sigma2_mu,
        common_variables, aux_identify_factor_variables
      )
      
      # Synchronize fitted values
      current_fit = get_predictions(curr_trees[j], x, single_tree = TRUE)
      yhat_bart = yhat_bart - tree_fits_store[, j] + current_fit
      tree_fits_store[, j] = current_fit
    }
    
    # 2. Update Error Model Parameters (Skew-t)
    y_hat     = yhat_bart
    residuals = y_scale - y_hat
    
    # Call the update controller for v, sigma2, lambda, and gamma2
    error_params = skewt_parameter_update(
      residuals, gamma2 = gamma2, v = v, d = d, 
      alpha_gamma = 2, beta_gamma = 1, 
      lambda = lambda, sigma_alpha = sigma_alpha, sigma_beta = sigma_beta
    )
    
    v      = error_params$v
    sigma2 = error_params$sigma2
    lambda = error_params$lambda
    gamma2 = error_params$gamma2
    
    # 3. Update Variable Selection Probabilities (Sparsity)
    if(isTRUE(sparse) && i > floor(total_iterations * 0.1)){
      s = update_s(var_count, p, 1)
    }
    
    # 4. Save Posterior Samples
    if ((i > nburn) && ((i - nburn) %% nthin == 0)) {
      curr = (i - nburn) / nthin
      tree_store[[curr]]      = curr_trees
      sigma2_store[curr]      = sigma2
      y_hat_store[curr, ]     = y_hat
      bart_store[curr, ]      = yhat_bart
      var_count_store[curr, ] = var_count
      s_prob_store[curr, ]    = s
      lambda_store[curr, ]    = lambda
      v_store[curr]           = v
      gamma2_store[curr]      = gamma2
    }
  }
  
  cat('\n MCMC Sampling Complete. \n')
  
  # --- Format Results (Rescale to original units) ---
  results <- list(
    trees            = tree_store,
    sigma2           = sigma2_store * y_sd^2,
    y_hat            = y_hat_store * y_sd + y_mean,
    bart_hat         = bart_store * y_sd + y_mean,
    v                = v_store,
    gamma2           = gamma2_store,
    npost            = npost,
    nburn            = nburn,
    nthin            = nthin,
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





