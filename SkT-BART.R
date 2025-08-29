#' @param y x
#' @param x x
#' @param sparse x
#' @param ntrees x
#' @param node_min_size x
#' @param alpha x
#' @param beta x
#' @param nu x
#' @param lambda x
#' @param mu_mu x
#' @param sigma2 x
#' @param sigma2_mu x
#' @param nburn x
#' @param npost x
#' @param nthin x
#'
#' @return x
#' @importFrom stats 'rgamma' 'rexp' 'dnorm' 'sd' 'rchisq' 'rnorm' 'pnorm' 'as.formula' 'terms' 'xtabs'
#' @importFrom truncnorm 'rtruncnorm'
#' @importFrom lme4 'lFormula'
#' @importFrom dbarts 'makeModelMatrixFromDataFrame'
#' @export
#'
#' @examples
#' #
#'
# sparse = FALSE
# ntrees = 50
# node_min_size = 5
# alpha = 0.95
# beta = 2
# nu = 3
# lambda = 0.1
# mu_mu = 0
# sigma2 = 1
# sigma2_mu = 1
# gamma2 = 1
# d = 0.1
# v = 1
# e = 0.5
# b = 0.5
# nburn = 100
# npost = 100
# nthin = 1
# lambda1 = runif(250, min = 0.1, max = 1)
bart_skewt = function(y,
                      x,
                      sparse = FALSE,
                      ntrees = 50,
                      node_min_size = 5,
                      gamma2 = 1,
                      alpha = 0.95,
                      beta = 2,
                      lambda1 = runif(250, min = 0.1, max = 1),
                      d = 0.1,
                      v = 1,
                      mu_mu = 0,
                      sigma2_mu = 1,
                      sigma2 = 1,
                      e = 0.5,
                      b = 0.5,
                      nburn = 1000,
                      npost = 1000,
                      nthin = 1) {
  
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
  lambda1_store    = matrix(NA, ncol = n, nrow = store_size)
  v_store          = rep(NA, store_size)
  
  # ---------- Standardize response variable ----------
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean) / y_sd
  
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
      lambda1_store[curr, ] = lambda1
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
      l_old = tree_full_skewt(
        tree = curr_trees[[j]],
        R = current_partial_residuals,
        lambda1 = lambda1,
        sigma2 = sigma2,
        sigma2_mu = sigma2_mu,
        gamma2 = gamma2,
        common_variables,
        aux_identify_factor_variables
      ) + get_tree_prior(curr_trees[[j]], alpha = alpha, beta = beta, common_variables)
      
      l_new = tree_full_skewt(
        tree = new_trees[[j]],
        R = current_partial_residuals,
        lambda1 = lambda1,
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
      curr_trees[[j]] = simulate_mu_skew(
        tree = curr_trees[[j]],
        R = current_partial_residuals,
        lambda1 = lambda1,
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
    hat = skewt_pro2_single(residuals, gamma2 = gamma2, v = v, d = d, alpha = alpha, beta = beta, lambda1 = lambda1)
    v       = mean(hat$v)
    sigma2  = mean(hat$sigma2)
    lambda1 = hat$lambda1
    gamma2  = mean(hat$gamma2)
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
                  lambda1=lambda1_store
  )
  class(results) <- "sktbart"
  return(results)
} # End main function
