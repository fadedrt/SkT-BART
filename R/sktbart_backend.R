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
List leaf_laplace_and_target_cpp(
    NumericVector R,
    NumericVector lambda,
    IntegerVector leaf_ids,
    IntegerVector unique_leaf_ids,
    NumericVector eval_mu_by_leaf,
    LogicalVector invalid_by_leaf,
    double gamma2,
    double sigma2,
    double sigma2_mu,
    bool draw
) {
  int N = R.size();
  int L = unique_leaf_ids.size();
  NumericVector mu_values(L);
  double log_q = 0.0;
  double log_kernel = 0.0;
  double gamma2_safe = std::max(gamma2, 1e-10);
  double pi = 3.14159265358979323846;

  std::map<int, int> leaf_to_pos;
  for (int l = 0; l < L; ++l) leaf_to_pos[unique_leaf_ids[l]] = l;

  std::vector< std::vector<int> > obs_by_leaf(L);
  for (int i = 0; i < N; ++i) {
    std::map<int, int>::iterator it = leaf_to_pos.find(leaf_ids[i]);
    if (it != leaf_to_pos.end()) obs_by_leaf[it->second].push_back(i);
  }

  for (int l = 0; l < L; ++l) {
    if (invalid_by_leaf[l]) {
      mu_values[l] = 0.0;
      continue;
    }

    const std::vector<int>& obs = obs_by_leaf[l];
    double mu_hat = 0.0;
    double laplace_var = sigma2_mu;

    if (!obs.empty()) {
      double sum_lam = 0.0;
      double sum_lam_r = 0.0;
      for (std::size_t k = 0; k < obs.size(); ++k) {
        int idx = obs[k];
        if (R_finite(R[idx]) && R_finite(lambda[idx])) {
          sum_lam += lambda[idx];
          sum_lam_r += lambda[idx] * R[idx];
        }
      }

      if (sum_lam > 0.0) {
        mu_hat = sum_lam_r / (sum_lam + sigma2 / sigma2_mu);

        for (int iter = 0; iter < 10; ++iter) {
          double old = mu_hat;
          double num = 0.0;
          double den = 0.0;

          for (std::size_t k = 0; k < obs.size(); ++k) {
            int idx = obs[k];
            if (R_finite(R[idx]) && R_finite(lambda[idx])) {
              double c_i = (R[idx] > mu_hat) ? (1.0 / gamma2_safe) : gamma2_safe;
              double w_i = c_i * lambda[idx];
              num += w_i * R[idx];
              den += w_i;
            }
          }

          mu_hat = num / (den + sigma2 / sigma2_mu);
          if (!R_finite(mu_hat)) mu_hat = old;
          if (std::abs(mu_hat - old) < 1e-4) break;
        }

        double sum_w = 0.0;
        for (std::size_t k = 0; k < obs.size(); ++k) {
          int idx = obs[k];
          if (R_finite(R[idx]) && R_finite(lambda[idx])) {
            double c_i = (R[idx] > mu_hat) ? (1.0 / gamma2_safe) : gamma2_safe;
            sum_w += c_i * lambda[idx];
          }
        }
        laplace_var = sigma2 / (sum_w + sigma2 / sigma2_mu);
      }
    }

    double laplace_sd = std::sqrt(std::max(laplace_var, 1e-14));
    double mu = draw ? R::rnorm(mu_hat, laplace_sd) : eval_mu_by_leaf[l];
    if (!R_finite(mu)) mu = mu_hat;
    mu_values[l] = mu;

    double z = (mu - mu_hat) / laplace_sd;
    log_q += -0.5 * std::log(2.0 * pi) - std::log(laplace_sd) - 0.5 * z * z;
    log_kernel += R::dnorm(mu, 0.0, std::sqrt(sigma2_mu), true);

    for (std::size_t k = 0; k < obs.size(); ++k) {
      int idx = obs[k];
      if (R_finite(R[idx]) && R_finite(lambda[idx])) {
        double resid = R[idx] - mu;
        double c_i = (R[idx] >= mu) ? (1.0 / gamma2_safe) : gamma2_safe;
        log_kernel += -0.5 * lambda[idx] * c_i * resid * resid / sigma2;
      }
    }
  }

  return List::create(
    Named("mu_values") = mu_values,
    Named("log_q") = log_q,
    Named("log_kernel") = log_kernel
  );
}
')

laplace_target_bundle <- function(tree, R, lambda, gamma2, sigma2, sigma2_mu,
                                  alpha, beta, common_vars, aux_factor_var,
                                  draw, eval_tree = NULL) {
  leaf_ids <- sort(unique(tree$node_indices))
  if (is.null(eval_tree)) eval_tree <- tree
  
  cpp <- leaf_laplace_and_target_cpp(
    R = R,
    lambda = lambda,
    leaf_ids = as.integer(tree$node_indices),
    unique_leaf_ids = as.integer(leaf_ids),
    eval_mu_by_leaf = eval_mu_by_leaf(eval_tree, leaf_ids),
    invalid_by_leaf = leaf_invalid_flags(tree, leaf_ids, common_vars, aux_factor_var),
    gamma2 = gamma2,
    sigma2 = sigma2,
    sigma2_mu = sigma2_mu,
    draw = draw
  )
  
  out_tree <- tree
  out_tree$tree_matrix[, "mu"] <- NA
  out_tree$tree_matrix[leaf_ids, "mu"] <- cpp$mu_values
  
  list(
    tree = out_tree,
    log_q = cpp$log_q,
    log_target = cpp$log_kernel + get_tree_prior(out_tree, alpha, beta, common_vars)
  )
}



