# ==============================================================================
# C++ Acceleration Backend for Skew-t BART (SkTBART)
# ==============================================================================
# Description:
# This script defines the high-performance C++ backend for the SkTBART model
# utilizing the Rcpp framework. It offloads computationally intensive bottleneck 
# tasks—such as tree topology traversals, data routing, and Iteratively Reweighted 
# Least Squares (IRLS) loops for Laplace approximations—from R to compiled C++.
#
# Performance Impact:
# By bypassing R's memory allocation and loop overheads, this module dramatically 
# accelerates MCMC sampling. This optimization is critical for scaling the algorithm 
# to handle high-dimensional, large-scale industrial datasets (e.g., semiconductor 
# manufacturing process control, chemical-mechanical planarization (CMP) polishing 
# liquid metrology, and real-time yield prediction).
# ==============================================================================

library(Rcpp)
# Compile C++ function directly within the R script
cppFunction('
List fill_tree_details_cpp(NumericMatrix tree_matrix, NumericMatrix X) {
    int n_obs = X.nrow();
    int n_nodes = tree_matrix.nrow();
    
    NumericMatrix new_tree_matrix = clone(tree_matrix);
    IntegerVector node_indices(n_obs, 1); 
    
    for(int i = 1; i < n_nodes; ++i) {
        int curr_parent = new_tree_matrix(i, 3) - 1; 
        int split_var = new_tree_matrix(curr_parent, 4) - 1; 
        double split_val = new_tree_matrix(curr_parent, 5);
        bool is_left = (new_tree_matrix(curr_parent, 1) - 1 == i);
        
        int node_size = 0;
        for(int obs = 0; obs < n_obs; ++obs) {
            if(node_indices[obs] == curr_parent + 1) {
                double x_val = X(obs, split_var);
                if(is_left) {
                    if(x_val < split_val) {
                        node_size++;
                        node_indices[obs] = i + 1;
                    }
                } else {
                    if(x_val >= split_val) {
                        node_size++;
                        node_indices[obs] = i + 1;
                    }
                }
            }
        }
        new_tree_matrix(i, 7) = node_size; 
    }
    
    return List::create(Named("tree_matrix") = new_tree_matrix,
                        Named("node_indices") = node_indices);
}
')

# Wrapper function to replace the original pure R loop
fill_tree_details = function(curr_tree, X) {
  # Call the compiled C++ backend
  res <- fill_tree_details_cpp(curr_tree$tree_matrix, X)
  
  return(list(tree_matrix = res$tree_matrix,
              node_indices = res$node_indices))
}



# Compile the C++ core acceleration code
cppFunction('
double laplace_irls_cpp(NumericVector R, 
                        NumericVector lambda, 
                        IntegerVector leaf_ids, 
                        IntegerVector unique_leaf_ids, 
                        NumericVector sigma2_mu_j, 
                        double gamma2, 
                        double sigma2) {
                        
    int N = R.size();
    int L = unique_leaf_ids.size();
    
    // Create a dictionary mapping Leaf ID to array index for fast accumulation
    std::map<int, int> leaf_to_idx;
    for(int l = 0; l < L; ++l) {
        leaf_to_idx[unique_leaf_ids[l]] = l;
    }
    
    NumericVector mu_hat(N, 0.0);
    double tol = 1e-4;
    int max_iter = 10;
    
    // Step 4: IRLS to find the local posterior mode
    for (int iter = 0; iter < max_iter; ++iter) {
        NumericVector mu_old = clone(mu_hat);
        NumericVector S_W(L, 0.0);
        NumericVector S_WR(L, 0.0);
        
        // Iterate over all observations, perform aggregate calculations
        for (int i = 0; i < N; ++i) {
            double c_i = (R[i] > mu_hat[i]) ? (1.0 / gamma2) : gamma2;
            double w_i = c_i * lambda[i];
            int idx = leaf_to_idx[leaf_ids[i]];
            S_W[idx] += w_i;
            S_WR[idx] += w_i * R[i];
        }
        
        double max_diff = 0.0;
        NumericVector mu_mode(L, 0.0);
        
        // Update mu for leaf nodes
        for (int l = 0; l < L; ++l) {
            if (sigma2_mu_j[l] == 0.0) {
                mu_mode[l] = 0.0;
            } else {
                double denom = S_W[l] + (sigma2 / sigma2_mu_j[l]);
                mu_mode[l] = S_WR[l] / denom;
            }
        }
        
        // Map mu back to the observation level and check for convergence
        for (int i = 0; i < N; ++i) {
            int idx = leaf_to_idx[leaf_ids[i]];
            mu_hat[i] = mu_mode[idx];
            double diff = std::abs(mu_hat[i] - mu_old[i]);
            if (diff > max_diff) max_diff = diff;
        }
        
        if (max_diff < tol) break;
    }
    
    // Step 5 and 6: Calculate the final statistics and log marginal likelihood at convergence
    NumericVector S_W(L, 0.0);
    NumericVector S_WR(L, 0.0);
    NumericVector S_WR2(L, 0.0);
    
    for (int i = 0; i < N; ++i) {
        double c_i = (R[i] > mu_hat[i]) ? (1.0 / gamma2) : gamma2;
        double w_i = c_i * lambda[i];
        int idx = leaf_to_idx[leaf_ids[i]];
        
        S_W[idx] += w_i;
        S_WR[idx] += w_i * R[i];
        S_WR2[idx] += w_i * R[i] * R[i];
    }
    
    double log_marg_lik = 0.0;
    for (int l = 0; l < L; ++l) {
        double denom = S_W[l] * sigma2_mu_j[l] + sigma2;
        double term1 = std::log(sigma2) - std::log(denom);
        double term2 = (sigma2_mu_j[l] * S_WR[l] * S_WR[l]) / (sigma2 * denom);
        double term3 = - S_WR2[l] / sigma2;
        log_marg_lik += 0.5 * (term1 + term2 + term3);
    }
    
    return log_marg_lik;
}
')

# A new, hybrid accelerated version of tree_full_skewt_laplace
tree_full_skewt_laplace <- function(tree, R, lambda, sigma2, sigma2_mu, gamma2, 
                                    common_vars, aux_factor_var) {
  
  # === Parts kept in R (handling complex topology validation) ===
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
  
  # === Parts delegated to C++ (handling high-density math operations) ===
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


library(Rcpp)

# Compile C++ core sampling acceleration code
cppFunction('
NumericVector simulate_mu_irls_cpp(NumericVector R, 
                                   NumericVector lambda, 
                                   IntegerVector leaf_ids, 
                                   IntegerVector unique_leaf_ids, 
                                   NumericVector sigma2_mu_j, 
                                   double gamma2, 
                                   double sigma2) {
    int N = R.size();
    int L = unique_leaf_ids.size();
    NumericVector mu_values(L, 0.0);
    
    // Create a dictionary: Leaf node ID -> array index
    std::map<int, int> leaf_to_idx;
    for(int l = 0; l < L; ++l) {
        leaf_to_idx[unique_leaf_ids[l]] = l;
    }
    
    // Highly efficient grouping: store the observation indices corresponding to each leaf node
    std::vector<std::vector<int>> obs_by_leaf(L);
    for(int i = 0; i < N; ++i) {
        int id = leaf_ids[i];
        if(leaf_to_idx.find(id) != leaf_to_idx.end()) {
            obs_by_leaf[leaf_to_idx[id]].push_back(i);
        }
    }
    
    double gamma2_safe = std::max(gamma2, 1e-10);
    double tol = 1e-4;
    
    // Perform independent calculations for each leaf node
    for(int l = 0; l < L; ++l) {
        double sig2_mu = sigma2_mu_j[l];
        if (sig2_mu == 0.0) {
            mu_values[l] = 0.0;
            continue;
        }
        
        const std::vector<int>& obs = obs_by_leaf[l];
        int n_leaf = obs.size();
        
        // If the leaf node has no valid observations, sample directly from the prior
        if (n_leaf == 0) {
            mu_values[l] = R::rnorm(0.0, std::sqrt(sig2_mu));
            continue;
        }
        
        // Initialize mu_hat
        double sum_lam = 0.0;
        double sum_lam_r = 0.0;
        for(int idx : obs) {
            if(std::isfinite(R[idx]) && std::isfinite(lambda[idx])) {
                sum_lam += lambda[idx];
                sum_lam_r += lambda[idx] * R[idx];
            }
        }
        sum_lam = std::max(sum_lam, 1e-10);
        double mu_hat = sum_lam_r / (sum_lam + sigma2 / sig2_mu);
        
        // Core IRLS loop (maximum 10 iterations)
        for(int iter = 0; iter < 10; ++iter) {
            double mu_old = mu_hat;
            double num = 0.0;
            double den = 0.0;
            for(int idx : obs) {
                if(std::isfinite(R[idx]) && std::isfinite(lambda[idx])) {
                    double R_i = R[idx];
                    double c_i = (R_i > mu_hat) ? (1.0 / gamma2_safe) : gamma2_safe;
                    double w_i = c_i * lambda[idx];
                    num += w_i * R_i;
                    den += w_i;
                }
            }
            mu_hat = num / (den + sigma2 / sig2_mu);
            if(std::abs(mu_hat - mu_old) < tol) break;
        }
        
        // Defensive programming: If mu_hat is abnormal (non-convergent or infinite), fall back to the median (consistent with original R code)
        if (!std::isfinite(mu_hat)) {
            std::vector<double> valid_R;
            for(int idx : obs) {
                if(std::isfinite(R[idx]) && std::isfinite(lambda[idx])) {
                    valid_R.push_back(R[idx]);
                }
            }
            if(!valid_R.empty()) {
                size_t m = valid_R.size() / 2;
                std::nth_element(valid_R.begin(), valid_R.begin() + m, valid_R.end());
                mu_hat = valid_R[m];
            } else {
                mu_hat = 0.0;
            }
        }
        
        // Calculate the Laplace variance and generate the final random sample
        double sum_w_final = 0.0;
        for(int idx : obs) {
            if(std::isfinite(R[idx]) && std::isfinite(lambda[idx])) {
                double R_i = R[idx];
                double c_i = (R_i > mu_hat) ? (1.0 / gamma2_safe) : gamma2_safe;
                sum_w_final += c_i * lambda[idx];
            }
        }
        
        double laplace_var = sigma2 / (sum_w_final + sigma2 / sig2_mu);
        double laplace_sd = std::sqrt(std::max(laplace_var, 1e-14));
        
        // R::rnorm seamlessly integrates with R's set.seed() system
        double mu_star = R::rnorm(mu_hat, laplace_sd);
        mu_values[l] = std::isfinite(mu_star) ? mu_star : mu_hat;
    }
    
    return mu_values;
}
')


# Re-wrapped pure R interface function
simulate_mu_skew_laplace <- function(tree, R, lambda1, gamma2, sigma2, sigma2_mu, 
                                     common_vars, aux_factor_var) {
  
  # 1. Handle complex tree node topology and multiple split validation logic at the R level (kept as is)
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
  
  # Align parameter dimensions: assign prior variance to all existing leaf nodes
  sigma2_mu_aux <- rep(sigma2_mu, length(unique_leaf_indices))
  
  # Forcibly set the variance of invalid nodes to 0
  if (length(which_terminal_no_double_split) > 0 && any(which_terminal_no_double_split != 0)) {
    invalid_idx <- which(unique_leaf_indices %in% which_terminal_no_double_split)
    sigma2_mu_aux[invalid_idx] <- 0
  }
  
  # 2. Call the C++ backend: strictly compress the time of the IRLS sampling loop
  mu_values <- simulate_mu_irls_cpp(
    R = R,
    lambda = lambda1,
    leaf_ids = as.integer(tree$node_indices),
    unique_leaf_ids = as.integer(unique_leaf_indices),
    sigma2_mu_j = sigma2_mu_aux,
    gamma2 = gamma2,
    sigma2 = sigma2
  )
  
  # 3. Write the calculated parameters back into the tree matrix intact
  tree$tree_matrix[, 'mu'] <- NA
  tree$tree_matrix[unique_leaf_indices, 'mu'] <- mu_values
  
  # Defensive overwrite, ensuring invalid nodes are forced to zero
  if (any(which_terminal_no_double_split != 0)) {
    tree$tree_matrix[which_terminal_no_double_split, 'mu'] <- 0 
  }
  
  return(tree)
}



