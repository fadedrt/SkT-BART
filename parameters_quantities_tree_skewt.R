get_tree_prior = function(tree, alpha, beta, common_vars) {

  # First find the level of each node, then the depth is the maximum of the level
  level = rep(NA, nrow(tree$tree_matrix))
  level[1] = 0 # First row always level 0

  # Escape quickly if tree is just a stump
  if(nrow(tree$tree_matrix) == 1) {
    return(log(1 - alpha)) # Tree depth is 0
  }

  for(i in 2:nrow(tree$tree_matrix)) {
    # Find the current parent
    curr_parent = as.numeric(tree$tree_matrix[i,'parent'])
    # This child must have a level one greater than it's current parent
    level[i] = level[curr_parent] + 1
  }

  # Only compute for the internal nodes
  internal_nodes = which(as.numeric(tree$tree_matrix[,'terminal']) == 0)
  log_prior = 0
  for(i in seq_along(internal_nodes)) {
    log_prior = log_prior + log(alpha) - beta * log1p(level[internal_nodes[i]])
  }
  # Now add on terminal nodes
  terminal_nodes = which(as.numeric(tree$tree_matrix[,'terminal']) == 1)
  for(i in seq_along(terminal_nodes)) {
    log_prior = log_prior + log(1 - alpha * ((1 + level[terminal_nodes[i]])^(-beta)))
  }

  return(log_prior)

  }

skewt_pro2_single <- function(x, gamma2, v, d = 0.15, alpha = 2, beta = 1, lambda) {
  n = length(x)
  
  # Indicator function: adjust scale based on sign
  indicator = ifelse(x >= 0, 1 / gamma2, gamma2)
  
  # Compute S for sigma2 update
  S = sum(lambda * x^2 * indicator)
  
  sigma2 = update_sigma2(S = S, n = n)
  v      = update_v_skewt(n = n, d = d, lambda = lambda)
  lambda = update_lambda_skewt(n, x, v, gamma2, sigma2 = sigma2)
  gamma2  = update_gamma2_skewt(n, lambda, x, sigma2 = sigma2, a = alpha, b = beta)
  
  return(list(
    v = v,
    gamma2 = gamma2,
    lambda = lambda,
    sigma2 = sigma2
  ))
}

update_gamma2_skewt <- function(n, lambda, residuals, sigma2, a, b) {
  # Log-target density function for sampling gamma2
  log_density <- function(gamma2) {
    gamma2 <- abs(gamma2)
    indicate <- ifelse(residuals >= 0, 1 / gamma2, gamma2)
    S <- sum(lambda * residuals^2 * indicate / sigma2)
    log_prior <- (n / 2 + a - 1) * log(gamma2) - n * log(gamma2 + 1)
    log_likelihood <- -S / 2 - b * gamma2
    return(log_prior + log_likelihood)
  }
  
  # Estimate lambda parameter for logf
  lambda1 <- find_lambda_one_d(log_density)
  
  # Bayesian sampling using ru package
  r <- ru(
    logf = log_density,
    d = 1,
    n = 50,
    trans = "BC",
    lambda = lambda1,
    lower = 0
  )
  
  return(mean(r$sim_vals))  # Return mean as the updated value
}


update_lambda_skewt <- function(n, residuals, v, gamma2, sigma2) {
  lambda <- numeric(n)
  for (i in seq_len(n)) {
    res_sq <- residuals[i]^2
    if (residuals[i] >= 0) {
      rate <- v / 2 + res_sq / (2 * gamma2 * sigma2)
    } else {
      rate <- v / 2 + res_sq * gamma2 / (2 * sigma2)
    }
    lambda[i] <- rgamma(1, shape = (v + 1) / 2, rate = rate)
  }
  return(lambda)
}

update_v_skewt <- function(n, d, lambda) {
  eta = sum(lambda - log(lambda)) / 2 + d
  T = n
  
  # Define the equation to solve for x*
  f <- function(x) {
    (T / 2) * (log(x / 2) + 1 - digamma(x / 2)) + (1 / x) - eta
  }
  
  # Solve for x* using uniroot
  solution <- uniroot(f, interval = c(1e-10, 100))
  x_star <- solution$root
  alpha <- 1 / x_star  # Set parameter for proposal distribution
  
  # Define kernel of the target distribution
  target_density <- function(x) {
    (x/2)^(T*x/2) * exp(-T * lgamma(x/2)) * exp(-eta * x)
  }
  
  # Define the proposal distribution
  proposal_density <- function(x) {
    alpha * exp(-alpha * x)
  }
  
  # Define acceptance probability function
  Q <- function(x) {
    exp((T*x/2) * log(x/2) - T * lgamma(x/2) - 
        (T*x_star/2) * log(x_star/2) + T * lgamma(x_star/2) + 
        (x_star - x) * eta - 1 + x / x_star)
  }
  
  # Accept-Reject sampling function
  ars_sample <- function(n) {
    samples <- numeric(n)
    accepted <- 0
    
    while (accepted < n) {
      # Sample from the proposal distribution
      x <- rexp(1, rate = alpha)
      
      # Compute acceptance probability
      prob <- Q(x)
      
      # Accept-Reject decision
      if (runif(1) < prob) {
        accepted <- accepted + 1
        samples[accepted] <- x
      }
    }
    
    return(samples)
  }
  
  return(ars_sample(1))
}

update_sigma2 <- function(S, n) {
  return(1 / rgamma(1, shape = n / 2, rate = S / 2))
}

tree_full_skewt <- function(tree, R, lambda, sigma2, sigma2_mu, gamma2 = gamma2, 
                            common_vars, aux_factor_var) {
  # Compute the log full conditional probability for the current tree (for MCMC acceptance)
  # tree: a single tree structure
  # R: current residual vector
  # lambda: latent weights for each observation (gamma samples)
  # sigma2: error variance
  # sigma2_mu: prior variance for mu in the tree
  # gamma2: skewness parameter
  # common_vars: set of covariates
  # aux_factor_var: set of factor variables that cannot be double split
  
  # -------------------------
  # Step 1: Identify all terminal nodes
  terminal_nodes <- which(tree$tree_matrix[,'terminal'] == 1)
  
  # -------------------------
  # Step 2: Check which terminals do not have double splits (or are stumps)
  if (nrow(tree$tree_matrix) != 1) {
    terminal_ancestors <- get_ancestors(tree)  # Get ancestors for all terminals
    ancestor_table <- table(terminal_ancestors[,1], terminal_ancestors[,2])
    
    terminals_invalid <- NULL
    for (k in 1:nrow(ancestor_table)) {
      ancestor_vars <- names(ancestor_table[k,])[ancestor_table[k, ] != 0]
      is_invalid <- any(sapply(aux_factor_var, function(x) all(ancestor_vars %in% x)))
      if (is_invalid) {
        terminals_invalid[k] <- rownames(ancestor_table)[k]
      }
    }
    
    # Keep invalid terminal node indices
    terminals_invalid <- as.numeric(terminals_invalid[!is.na(terminals_invalid)])
  } else {
    terminals_invalid <- 0  # stump case
  }
  
  # -------------------------
  # Step 3: Set sigma2_mu = 0 for invalid terminal nodes
  sigma2_mu_vec <- rep(sigma2_mu, length(terminal_nodes))
  sigma2_mu_vec[which(terminal_nodes %in% terminals_invalid)] <- 0
  
  # -------------------------
  # Step 4: Compute information for each terminal node
  node_sizes <- as.numeric(tree$tree_matrix[terminal_nodes, 'node_size'])
  
  # Compute asymmetric skewness weight factor C_i
  C_i <- 1 / gamma2 * pnorm(R / sigma2_mu) + gamma2 * (1 - pnorm(R / sigma2_mu))
  
  # Sum over each node
  S0 <- C_i * lambda * R^2
  S_j2 <- tapply(S0, tree$node_indices, sum)
  
  S1 <- lambda * R * C_i
  S_j1 <- tapply(S1, tree$node_indices, sum)
  
  M_j <- tapply(lambda, tree$node_indices, sum)
  
  # -------------------------
  # Step 5: Compute log posterior probability
  denom <- M_j * sigma2_mu_vec + sigma2
  log_post <- 0.5 * (
    sum(log(sigma2) - log(denom)) +
      sum((sigma2_mu_vec * S_j1^2) / (sigma2 * denom))
  )
  
  return(log_post)
}

draw_mu_by_mh <- function(R, lambda, gamma2, sigma2, sigma2_mu,
                          init_mu = median(R),
                          proposal_sd = mad(R) / 2,
                          max_iter = 500) {
  # Metropolis-Hastings sampling for a single mu (used in Skew-t BART)
  
  log_target <- function(mu) {
    c <- ifelse(R - mu >= 0, 1 / gamma2, gamma2)
    -sum(lambda * c * (R - mu)^2) / (2 * sigma2) - mu^2 / (2 * sigma2_mu)
  }
  
  mu <- init_mu
  samples <- numeric(max_iter)
  acc <- 0
  
  for (i in seq_len(max_iter)) {
    prop <- rnorm(1, mu, proposal_sd)
    log_alpha <- log_target(prop) - log_target(mu)
    if (log(runif(1)) < log_alpha) {
      mu <- prop
      acc <- acc + 1
    }
    samples[i] <- mu
    
    # Optional: adapt proposal_sd to keep acceptance rate between 20%-60%
    if (i %% 50 == 0 && i >= 100) {
      rate <- acc / i
      if (rate >= 0.2 && rate <= 0.6) break
      proposal_sd <- proposal_sd * if (rate < 0.2) 1.5 else 0.7
    }
  }
  
  return(mu)
}


simulate_mu_skewt = function(tree, R, lambda, gamma2, sigma2, sigma2_mu, 
                            common_vars, aux_factor_var) {
  # Simulate mu values for a given tree
  # First find which rows are terminal nodes
  which_terminal = which(tree$tree_matrix[, 'terminal'] == 1)
  # Identify those terminals that don't have a double split on g and e
  if (nrow(tree$tree_matrix) != 1) {
    terminal_ancestors = get_ancestors(tree)  # get the ancestor for all terminal nodes
    aux_table = table(terminal_ancestors[, 1], terminal_ancestors[, 2])  # create a table
    which_terminal_no_double_split = NULL
    for (k in 1:nrow(aux_table)) {
      split_var_ancestors = names(aux_table[k, ])[aux_table[k, ] != 0]
      invalid_ancestors = any(unlist(
        lapply(aux_factor_var, function(x) all(split_var_ancestors %in% x))
      ))
      if (invalid_ancestors == TRUE) {
        which_terminal_no_double_split[k] = rownames(aux_table)[k]
      }
    }
    # terminals with only one ancestor where the ancestor is common to X1 and X2
    which_terminal_no_double_split = as.numeric(
      which_terminal_no_double_split[!is.na(which_terminal_no_double_split)]
    )
  } else {
    which_terminal_no_double_split = 0  # stump
  }
  # set up sigma2_mu = 0 for all which_terminal_no_double_split
  sigma2_mu_aux = rep(sigma2_mu, length(which_terminal))
  sigma2_mu_aux[which(which_terminal %in% which_terminal_no_double_split)] = 0
  nj = as.numeric(tree$tree_matrix[which_terminal, 'node_size'])
  mu <- c()
  for (i in 1:length(nj)) {
    ind <- which(tree$node_indices == sort(unique(tree$node_indices))[i])
    mu[i] <- draw_mu_by_mh(
      R = R[ind], lambda = lambda[ind], gamma2 = gamma2, sigma2 = sigma2, 
      sigma2_mu = sigma2_mu, init_mu = median(R), proposal_sd = mad(R) / 2
    )
  }
  # Wipe all the old mus out for other nodes
  tree$tree_matrix[, 'mu'] = NA
  # Put in just the ones that are useful
  tree$tree_matrix[which_terminal, 'mu'] = mu
  tree$tree_matrix[which_terminal_no_double_split, 'mu'] = 0  # set to zero the terminal node with no interaction
  return(tree)
}







