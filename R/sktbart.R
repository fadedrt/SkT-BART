#' BART with Skew-t Errors via Laplace proposal
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
#' @importFrom skewt rskt qskt
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
  if (is.null(sigma2_mu)) sigma2_mu <- 0.5^2 / ntrees
  if (is.null(lambda)) lambda <- rep(1, length(y))
  
  aux_identify_factor_variables <- NULL
  common_variables <- NULL
  
  total_iterations <- nburn + npost * nthin
  n <- length(y)
  p <- ncol(x)
  
  tree_store <- vector("list", npost)
  gamma2_store <- rep(NA, npost)
  sigma2_store <- rep(NA, npost)
  y_hat_store <- matrix(NA, ncol = n, nrow = npost)
  bart_store <- matrix(NA, ncol = n, nrow = npost)
  var_count <- rep(0, p)
  var_count_store <- matrix(0, ncol = p, nrow = npost)
  s_prob_store <- matrix(0, ncol = p, nrow = npost)
  v_store <- rep(NA, npost)
  lambda_store <- matrix(NA, ncol = n, nrow = npost)
  tree_accept_store <- rep(NA, npost)
  
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
    accepted_this_iter <- 0
    
    for (j in seq_len(ntrees)) {
      current_partial_residuals <- y_scale - yhat_bart + tree_fits_cache[, j]
      
      upd <- joint_update_tree_leaf_fast(
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
        move_type <- upd$move_type
        if (move_type == "change" && length(curr_trees[[j]]$var) >= 2) {
          new_var <- curr_trees[[j]]$var[1]
          old_var <- curr_trees[[j]]$var[2]
          if (new_var %in% seq_len(p)) var_count[new_var] <- var_count[new_var] + 1
          if (old_var %in% seq_len(p)) var_count[old_var] <- var_count[old_var] - 1
        }
        if (move_type == "grow" && length(curr_trees[[j]]$var) > 0) {
          for (vv in curr_trees[[j]]$var) {
            if (vv %in% seq_len(p)) var_count[vv] <- var_count[vv] + 1
          }
        }
        if (move_type == "prune" && length(curr_trees[[j]]$var) > 0) {
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
  cat("\\n MCMC Sampling Complete. \\n")
  
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


