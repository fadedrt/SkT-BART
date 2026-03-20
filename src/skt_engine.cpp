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
    
    // 建立 叶子节点ID 到 数组索引 的映射字典，方便快速累加
    std::map<int, int> leaf_to_idx;
    for(int l = 0; l < L; ++l) {
        leaf_to_idx[unique_leaf_ids[l]] = l;
    }
    
    NumericVector mu_hat(N, 0.0);
    double tol = 1e-4;
    int max_iter = 10;
    
    // 第4步：IRLS 寻找局部后验众数
    for (int iter = 0; iter < max_iter; ++iter) {
        NumericVector mu_old = clone(mu_hat);
        NumericVector S_W(L, 0.0);
        NumericVector S_WR(L, 0.0);
        
        // 遍历所有观测值，聚合计算
        for (int i = 0; i < N; ++i) {
            double c_i = (R[i] > mu_hat[i]) ? (1.0 / gamma2) : gamma2;
            double w_i = c_i * lambda[i];
            int idx = leaf_to_idx[leaf_ids[i]];
            S_W[idx] += w_i;
            S_WR[idx] += w_i * R[i];
        }
        
        double max_diff = 0.0;
        NumericVector mu_mode(L, 0.0);
        
        // 更新叶子节点的 mu
        for (int l = 0; l < L; ++l) {
            if (sigma2_mu_j[l] == 0.0) {
                mu_mode[l] = 0.0;
            } else {
                double denom = S_W[l] + (sigma2 / sigma2_mu_j[l]);
                mu_mode[l] = S_WR[l] / denom;
            }
        }
        
        // 将 mu 映射回观测值层面并检查收敛
        for (int i = 0; i < N; ++i) {
            int idx = leaf_to_idx[leaf_ids[i]];
            mu_hat[i] = mu_mode[idx];
            double diff = std::abs(mu_hat[i] - mu_old[i]);
            if (diff > max_diff) max_diff = diff;
        }
        
        if (max_diff < tol) break;
    }
    
    // 第5步和第6步：在收敛状态下计算最终统计量与对数边缘似然
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
    
    // 建立字典：叶子节点 ID -> 数组索引
    std::map<int, int> leaf_to_idx;
    for(int l = 0; l < L; ++l) {
        leaf_to_idx[unique_leaf_ids[l]] = l;
    }
    
    // 极其高效的分组：把每个叶子节点对应的观测值索引存起来
    std::vector<std::vector<int>> obs_by_leaf(L);
    for(int i = 0; i < N; ++i) {
        int id = leaf_ids[i];
        if(leaf_to_idx.find(id) != leaf_to_idx.end()) {
            obs_by_leaf[leaf_to_idx[id]].push_back(i);
        }
    }
    
    double gamma2_safe = std::max(gamma2, 1e-10);
    double tol = 1e-4;
    
    // 对每个叶子节点进行独立计算
    for(int l = 0; l < L; ++l) {
        double sig2_mu = sigma2_mu_j[l];
        if (sig2_mu == 0.0) {
            mu_values[l] = 0.0;
            continue;
        }
        
        const std::vector<int>& obs = obs_by_leaf[l];
        int n_leaf = obs.size();
        
        // 如果该叶节点没有有效观测值，直接从先验抽样
        if (n_leaf == 0) {
            mu_values[l] = R::rnorm(0.0, std::sqrt(sig2_mu));
            continue;
        }
        
        // 初始化 mu_hat
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
        
        // IRLS 核心循环 (最多迭代10次)
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
        
        // 防御性编程：如果 mu_hat 出现异常（不收敛或无穷大），回退到中位数（和原 R 代码保持一致）
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
        
        // 计算 Laplace 方差并生成最终的随机样本
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
        
        // R::rnorm 可以无缝对接 R 语言的 set.seed() 系统
        double mu_star = R::rnorm(mu_hat, laplace_sd);
        mu_values[l] = std::isfinite(mu_star) ? mu_star : mu_hat;
    }
    
    return mu_values;
}
')

