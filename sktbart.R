#' BART with Skew-t Errors via Laplace Approximation
#'
#' @param y A numeric vector of the response variable.
#' @param x A data frame or matrix of predictors.
#' @param sparse Logical; if TRUE, uses a Dirichlet prior for sparse variable selection.
#' @param ntrees Number of trees in the ensemble (default is 50).
#' @param node_min_size Minimum number of observations required in a terminal node.
#' @param alpha Base hyperparameter for tree prior (controls depth).
#' @param beta Power hyperparameter for tree prior (controls depth).
#' @param v Degrees of freedom for the Skew-t distribution (initial value).
#' @param lambda Latent weights for the scale-mixture representation (initial vector).
#' @param mu_mu Mean of the prior for leaf node parameters.
#' @param sigma2 Initial error variance.
#' @param sigma2_mu Prior variance for leaf node parameters.
#' @param nburn Number of burn-in iterations for MCMC.
#' @param npost Number of posterior samples to collect.
#' @param nthin Thinning interval for MCMC.
#'
#' @return An object of class 'sktbart' containing posterior samples.
#' 
#' @importFrom stats rgamma rexp dnorm sd rchisq rnorm pnorm as.formula terms xtabs lm var median
#' @importFrom truncnorm rtruncnorm
#' @importFrom dbarts makeModelMatrixFromDataFrame
#' @importFrom rust ru
#' @export
#'
#' @examples
#' # Example: Friedman Nonlinear Model
#' n <- 250; p <- 10
#' x <- matrix(runif(n * p), n, p)
#' # Formula: 10*sin(pi*x1*x2) + 20*(x3-0.5)^2 + 10*x4 + 5*x5
#' f_x <- 10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 10 * x[,4] + 5 * x[,5]
#' y <- f_x + rnorm(n, sd = 1.0)
#' 
#' model <- sktbart(y, x, nburn = 2500, npost = 2500)
#'
bart_skewt_laplace <- function(y,
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
                       sigma_alpha=1.5
) {
  if (is.null(sigma2_mu)) {
    sigma2_mu = 0.5^2 / ntrees  
  }
  if (is.null(lambda)){
    lambda =rep(1,length(y))
  }
  aux_identify_factor_variables = NULL
  common_variables = NULL
  
  TotIter = nburn + npost * nthin
  store_size = npost
  n = length(y)
  p = ncol(x)
  
  # Initialize storage
  tree_store       = vector('list', store_size)
  gamma2_store     = rep(NA, store_size)
  sigma2_store     = rep(NA, store_size)
  y_hat_store      = matrix(NA, ncol = n, nrow = store_size)
  bart_store       = matrix(NA, ncol = n, nrow = store_size)
  var_count        = rep(0, p)
  var_count_store  = matrix(0, ncol = p, nrow = store_size)
  s_prob_store     = matrix(0, ncol = p, nrow = store_size)
  tree_fits_store  = matrix(0, ncol = ntrees, nrow = n)
  lambda_store    = matrix(NA, ncol = n, nrow = store_size)
  v_store          = rep(NA, store_size)
  
  # ---------- Standardize response variable ----------
  y_mean = mean(y)
  y_sd = sd(y) 
  y_scale = (y - y_mean) / y_sd
  lm_fit = lm(y_scale ~ x)
  sigma_hat_sq = var(lm_fit$residuals)
  sigma_beta  = 0.292 * sigma_hat_sq 
  
  # Initialize variable inclusion probabilities
  s = rep(1 / p, p)
  current_partial_residuals = y_scale
  
  # Initialize tree stumps
  curr_trees = create_stump(num_trees = ntrees, y = y_scale, X = x)
  new_trees  = curr_trees
  yhat_bart  = get_predictions(curr_trees, x, single_tree = (ntrees == 1))
  
  pb = utils::txtProgressBar(min = 1, max = TotIter, style = 3, width = 60)
  
  for (i in seq_len(TotIter)) {
    utils::setTxtProgressBar(pb, i)
    
    # Store posterior draws after burn-in and thinning
    if ((i > nburn) && ((i - nburn) %% nthin == 0)) {
      curr = (i - nburn) / nthin
      tree_store[[curr]]    = curr_trees
      sigma2_store[curr]    = sigma2
      y_hat_store[curr, ]   = y_hat
      bart_store[curr, ]    = yhat_bart
      var_count_store[curr, ] = var_count
      s_prob_store[curr, ]  = s
      lambda_store[curr, ] = lambda
      v_store[curr]         = v
      gamma2_store[curr]    = gamma2
    }
    
    for (j in seq_len(ntrees)) {
      current_partial_residuals = y_scale - yhat_bart + tree_fits_store[, j]
      
      type = sample_move(curr_trees[[j]], i, nburn)
      new_trees[[j]] = update_tree(
        y = current_partial_residuals, X = x,
        type = type, curr_tree = curr_trees[[j]],
        node_min_size = node_min_size, s = s,
        common_vars = common_variables,
        aux_factor_var = aux_identify_factor_variables
      )
      
      # Compute full tree likelihood
      l_old = tree_full_skewt_laplace(
        tree = curr_trees[[j]],
        R = current_partial_residuals,
        lambda = lambda,
        sigma2 = sigma2,
        sigma2_mu = sigma2_mu,
        gamma2 = gamma2,
        common_variables,
        aux_identify_factor_variables
      ) + get_tree_prior(curr_trees[[j]], alpha = alpha, beta = beta, common_variables)
      
      l_new = tree_full_skewt_laplace(
        tree = new_trees[[j]],
        R = current_partial_residuals,
        lambda = lambda,
        sigma2 = sigma2,
        sigma2_mu = sigma2_mu,
        gamma2 = gamma2,
        common_variables,
        aux_identify_factor_variables
      ) + get_tree_prior(new_trees[[j]], alpha = alpha, beta = beta, common_variables)
      
      a = if (isTRUE(new_trees[[j]]$ForceStump)) 1 else l_new - l_old
      if (a > 0 || a > log(runif(1))) {
        curr_trees[[j]] = new_trees[[j]]
        if (type == 'change') {
          var_count[curr_trees[[j]]$var[1]] = var_count[curr_trees[[j]]$var[1]] + 1
          var_count[curr_trees[[j]]$var[2]] = var_count[curr_trees[[j]]$var[2]] - 1
        }
        if (type == 'grow') var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] + 1
        if (type == 'prune') var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] - 1
      }
      
      # Update terminal node parameters
      curr_trees[[j]] = simulate_mu_skew_laplace(
        tree = curr_trees[[j]],
        R = current_partial_residuals,
        lambda = lambda,
        gamma2 = gamma2,
        sigma2 = sigma2,
        sigma2_mu = sigma2_mu,
        common_variables,
        aux_identify_factor_variables
      )
      
      # Update fitted values
      current_fit = get_predictions(curr_trees[j], x, single_tree = TRUE)
      yhat_bart = yhat_bart - tree_fits_store[, j] + current_fit
      tree_fits_store[, j] = current_fit
    }
    
    y_hat = yhat_bart
    residuals = y_scale - y_hat
    
    # Update skew-t parameters
    hat = skewt_pro2_single_prior(residuals, gamma2 = gamma2, v = v, d = d, alpha_gamma = 2, beta_gamma = 1, lambda = lambda, sigma_alpha, sigma_beta)#
    v       = hat$v
    sigma2  = hat$sigma2
    lambda = hat$lambda
    gamma2  = hat$gamma2
    
    if(isTRUE(sparse) & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }
  }
  
  
  cat('\n')
  results <- list(trees = tree_store,
                  sigma2 = sigma2_store*y_sd^2,
                  y_hat = y_hat_store*y_sd + y_mean,
                  bart_hat = bart_store*y_sd+ y_mean,
                  v=v_store,
                  gamma2=gamma2_store,
                  npost = npost,
                  nburn = nburn,
                  nthin = nthin,
                  ntrees = ntrees,
                  y_mean = y_mean,
                  y_sd = y_sd,
                  var_count_store = var_count_store,
                  s = s_prob_store,
                  lambda=lambda_store
  )
  class(results) <- "sktbart"
  return(results)
} # End main function




