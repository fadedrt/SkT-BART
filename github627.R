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
                      nburn = 2500,
                      npost = 2500,
                      nthin = 1,
                      robust_scaling = FALSE) {
  
  aux_identify_factor_variables = NULL
  common_variables = NULL
  
  TotIter = nburn + npost * nthin
  store_size = npost
  n = length(y)
  p = ncol(x)
  
  # 初始化存储器
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
  
  # ---------- 响应变量标准化 ----------
  if (robust_scaling) {
    y_mean = median(y)
    y_sd = IQR(y)
  } else {
    y_mean = mean(y)
    y_sd = sd(y)
  }
  y_scale = (y - y_mean) / y_sd
  
  # 初始化
  s = rep(1 / p, p)
  current_partial_residuals = y_scale
  
  # 初始化树桩
  curr_trees = create_stump(num_trees = ntrees, y = y_scale, X = x)
  new_trees  = curr_trees
  yhat_bart  = get_predictions(curr_trees, x, single_tree = (ntrees == 1))
  
  pb = utils::txtProgressBar(min = 1, max = TotIter, style = 3, width = 60)
  
  for (i in seq_len(TotIter)) {
    utils::setTxtProgressBar(pb, i)
    
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
      
      # 计算全树似然
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
      
      # 更新每棵树叶节点参数
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
      
      # 更新拟合值
      current_fit = get_predictions(curr_trees[j], x, single_tree = TRUE)
      yhat_bart = yhat_bart - tree_fits_store[, j] + current_fit
      tree_fits_store[, j] = current_fit
    }
    
    y_hat = yhat_bart
    residuals = y_scale - y_hat
    
    hat = skewt_pro2_single(residuals, gamma2 = gamma2, v = v, d = d, alpha = alpha, beta = beta, lambda1 = lambda1)
    v       = mean(hat$v)
    sigma2  = mean(hat$sigma2)
    lambda1 = hat$lambda1
    gamma2  = mean(hat$gamma2)
  }
  
  cat('\n')
  
  return(structure(list(
    trees = tree_store,
    sigma2 = sigma2_store * y_sd^2,
    y_hat = y_hat_store * y_sd + y_mean,
    bart_hat = bart_store * y_sd + y_mean,
    v = v_store,
    gamma2 = gamma2_store,
    npost = npost,
    nburn = nburn,
    nthin = nthin,
    ntrees = ntrees,
    y_mean = y_mean,
    y_sd = y_sd,
    var_count_store = var_count_store,
    s = s_prob_store,
    lambda1 = lambda1_store
  ), class = "bart"))
}



skewt_pro2_single <- function(x, gamma2, v, d = 0.15, alpha = 2, beta = 1, lambda1) {
  n = length(x)
  
  # 指示函数：根据正负调整尺度
  indicator = ifelse(x >= 0, 1 / gamma2, gamma2)
  
  # 计算 S 用于 sigma2 更新
  S = sum(lambda1 * x^2 * indicator)
  
  sigma2 = update_sigma2(S = S, n = n)
  v      = update_v_skewt(n = n, d = d, lambda1 = lambda1)
  lambda1 = update_lambda1_skewt(n, x, v, gamma2, sigma2 = sigma2)
  gamma2  = update_gamma2_skewt(n, lambda1, x, sigma2 = sigma2, a = alpha, b = beta)
  
  return(list(
    v = v,
    gamma2 = gamma2,
    lambda1 = lambda1,
    sigma2 = sigma2
  ))
}


update_gamma2_skewt <- function(n, lambda1, residuals, sigma2, a, b) {
  # 对数目标密度函数，用于 gamma2 的采样
  log_density <- function(gamma2) {
    gamma2 <- abs(gamma2)
    indicate <- ifelse(residuals >= 0, 1 / gamma2, gamma2)
    S <- sum(lambda1 * residuals^2 * indicate / sigma2)
    log_prior <- (n / 2 + a - 1) * log(gamma2) - n * log(gamma2 + 1)
    log_likelihood <- -S / 2 - b * gamma2
    return(log_prior + log_likelihood)
  }
  
  # 估计 logf 的 lambda 参数
  lambda <- find_lambda_one_d(log_density)
  
  # 使用 ru 包进行贝叶斯采样
  r <- ru(
    logf = log_density,
    d = 1,
    n = 50,
    trans = "BC",
    lambda = lambda,
    lower = 0
  )
  
  return(mean(r$sim_vals))  # 返回平均值作为更新值
}


update_lambda1_skewt <- function(n, residuals, v, gamma2, sigma2) {
  lambda1 <- numeric(n)
  for (i in seq_len(n)) {
    res_sq <- residuals[i]^2
    if (residuals[i] >= 0) {
      rate <- v / 2 + res_sq / (2 * gamma2 * sigma2)
    } else {
      rate <- v / 2 + res_sq * gamma2 / (2 * sigma2)
    }
    lambda1[i] <- rgamma(1, shape = (v + 1) / 2, rate = rate)
  }
  return(lambda1)
}


update_v_skewt <- function(n, d, lambda1) {
  eta <- sum(lambda1 - log(lambda1)) / 2 + d
  T <- n
  
  # 方程用于解出 x_star
  f <- function(x) {
    (T / 2) * (log(x / 2) + 1 - digamma(x / 2)) + (1 / x) - eta
  }
  
  # 求解 x_star
  solution <- uniroot(f, interval = c(1e-8, 100))
  x_star <- solution$root
  alpha <- 1 / x_star
  
  # 目标分布核函数（非标准化）
  target_log_density <- function(x) {
    (T * x / 2) * log(x / 2) - T * lgamma(x / 2) - eta * x
  }
  
  # 接受概率函数（与 proposal 比例）
  Q <- function(x) {
    exp(target_log_density(x) - target_log_density(x_star) + (x_star - x) * eta - 1 + x / x_star)
  }
  
  # 接受-拒绝采样
  ars_sample <- function(n) {
    samples <- numeric(n)
    accepted <- 0
    while (accepted < n) {
      x <- rexp(1, rate = alpha)
      if (runif(1) < Q(x)) {
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


tree_full_skewt <- function(tree, R, lambda1, sigma2, sigma2_mu, gamma2 = gamma2, 
                            common_vars, aux_factor_var) {
  # 计算当前树的 log full conditional 概率，用于 MCMC 样本接受判断
  # tree: 一棵单独的树结构
  # R: 当前残差向量
  # lambda1: 每个观测对应的 latent weight（gamma分布样本）
  # sigma2: 误差方差
  # sigma2_mu: 树结构中 mu 的先验方差
  # gamma2: 偏态参数
  # common_vars: 共变量变量集合
  # aux_factor_var: 不允许双重分裂的因子变量集合
  
  # -------------------------
  # Step 1: 找到所有 terminal nodes（叶子节点）
  terminal_nodes <- which(tree$tree_matrix[,'terminal'] == 1)
  
  # -------------------------
  # Step 2: 判断哪些 terminal 没有 double split（或者是 stump）
  if (nrow(tree$tree_matrix) != 1) {
    terminal_ancestors <- get_ancestors(tree)  # 获取所有 terminal 的祖先节点
    ancestor_table <- table(terminal_ancestors[,1], terminal_ancestors[,2])
    
    terminals_invalid <- NULL
    for (k in 1:nrow(ancestor_table)) {
      ancestor_vars <- names(ancestor_table[k,])[ancestor_table[k, ] != 0]
      is_invalid <- any(sapply(aux_factor_var, function(x) all(ancestor_vars %in% x)))
      if (is_invalid) {
        terminals_invalid[k] <- rownames(ancestor_table)[k]
      }
    }
    
    # 保留非法节点编号
    terminals_invalid <- as.numeric(terminals_invalid[!is.na(terminals_invalid)])
  } else {
    terminals_invalid <- 0  # stump 情况
  }
  
  # -------------------------
  # Step 3: 设置 sigma2_mu = 0 对于非法 terminal 节点
  sigma2_mu_vec <- rep(sigma2_mu, length(terminal_nodes))
  sigma2_mu_vec[which(terminal_nodes %in% terminals_invalid)] <- 0
  
  # -------------------------
  # Step 4: 计算每个 terminal 节点的信息量
  node_sizes <- as.numeric(tree$tree_matrix[terminal_nodes, 'node_size'])
  
  # 计算非对称偏态权重因子 C_i
  C_i <- 1 / gamma2 * pnorm(R / sigma2_mu) + gamma2 * (1 - pnorm(R / sigma2_mu))
  
  # 对每个节点分组求和
  S0 <- C_i * lambda1 * R^2
  S_j2 <- tapply(S0, tree$node_indices, sum)
  
  S1 <- lambda1 * R * C_i
  S_j1 <- tapply(S1, tree$node_indices, sum)
  
  M_j <- tapply(lambda1, tree$node_indices, sum)
  
  # -------------------------
  # Step 5: 计算对数后验概率
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
  # 用 Metropolis-Hastings 采样单个 mu 值（用于 Skew-t BART）
  
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
    
    # 可选：自动调整 proposal_sd，以控制接受率在 20%-60% 区间
    if (i %% 50 == 0 && i >= 100) {
      rate <- acc / i
      if (rate >= 0.2 && rate <= 0.6) break
      proposal_sd <- proposal_sd * if (rate < 0.2) 1.5 else 0.7
    }
  }
  
  return(mu)
}


simulate_mu_skew <- function(tree, R, lambda1, gamma2, sigma2, sigma2_mu,
                             common_vars, aux_factor_var) {
  # 为树中的 terminal node 根据偏态误差结构采样 mu 值
  
  terminal_nodes <- which(tree$tree_matrix[, 'terminal'] == 1)
  
  # 检查 terminal node 是否含有不允许的双重分裂路径
  if (nrow(tree$tree_matrix) != 1) {
    terminal_ancestors <- get_ancestors(tree)
    ancestor_table <- table(terminal_ancestors[, 1], terminal_ancestors[, 2])
    
    invalid_terminals <- NULL
    for (k in 1:nrow(ancestor_table)) {
      ancestor_vars <- names(ancestor_table[k, ])[ancestor_table[k, ] != 0]
      is_invalid <- any(sapply(aux_factor_var, function(x) all(ancestor_vars %in% x)))
      if (is_invalid) {
        invalid_terminals[k] <- rownames(ancestor_table)[k]
      }
    }
    
    invalid_terminals <- as.numeric(invalid_terminals[!is.na(invalid_terminals)])
  } else {
    invalid_terminals <- 0  # stump 情况
  }
  
  # 针对不允许交互的 terminal，设置 mu 的先验方差为 0
  sigma2_mu_vec <- rep(sigma2_mu, length(terminal_nodes))
  sigma2_mu_vec[terminal_nodes %in% invalid_terminals] <- 0
  
  node_sizes <- as.numeric(tree$tree_matrix[terminal_nodes, 'node_size'])
  
  # 对每个 terminal node 独立采样 mu
  mu <- numeric(length(node_sizes))
  unique_indices <- sort(unique(tree$node_indices))  # 避免多次排序
  
  for (i in seq_along(node_sizes)) {
    idx <- which(tree$node_indices == unique_indices[i])
    mu[i] <- draw_mu_by_mh(
      R = R[idx],
      lambda = lambda1[idx],
      gamma2 = gamma2,
      sigma2 = sigma2,
      sigma2_mu = sigma2_mu_vec[i],
      init_mu = median(R[idx]),
      proposal_sd = mad(R[idx]) / 2
    )
  }
  
  # 清除旧 mu，填入新的 mu 值
  tree$tree_matrix[, 'mu'] <- NA
  tree$tree_matrix[terminal_nodes, 'mu'] <- mu
  tree$tree_matrix[invalid_terminals, 'mu'] <- 0  # 直接来自根节点的“禁止交互”终端节点
  
  return(tree)
}


#### -----------------------    初始树       ----------------------------####
create_stump = function(num_trees,
                        y,
                        X) {
  
  # Each tree is a list of 2 elements
  # The 2 elements are the tree matrix (8 columns), and the node indices
  # The columns of the tree matrix are:
  # Terminal (0 = no, 1 = yes)
  # Child left
  # Child right
  # Node parents
  # Split variable
  # Split value
  # mu
  # Node size
  
  # Create holder for trees
  all_trees = vector('list', length = num_trees)
  # Loop through trees
  for (j in 1:num_trees) {
    # Set up each tree to have two elements in the list as described above
    all_trees[[j]] = vector('list', length = 2)
    # Give the elements names
    names(all_trees[[j]]) = c('tree_matrix',
                              'node_indices')
    # Create the two elements: first is a matrix
    all_trees[[j]][[1]] = matrix(NA, ncol = 8, nrow = 1)
    
    # Second is the assignment to node indices
    all_trees[[j]][[2]] = rep(1, length(y))
    
    # Create column names
    colnames(all_trees[[j]][[1]]) = c('terminal',
                                      'child_left',
                                      'child_right',
                                      'parent',
                                      'split_variable',
                                      'split_value',
                                      'mu',
                                      'node_size')
    
    # Set values for stump
    all_trees[[j]][[1]][1,] = c(1, NA, NA, NA, NA, NA, 0 , length(y))
    all_trees[[j]]$ForceStump = FALSE
    
  } # End of loop through trees
  
  return(all_trees)
  
}
#### -----------------------  树的更新        ----------------------------####
update_tree = function(y, # Target variable
                       X, # Feature matrix
                       type = c('grow',   # Grow existing tree
                                'prune',  # Prune existing tree
                                'change', # Change existing tree - change split variable and value for an internal node
                                'swap'),  # Swap existing tree - swap splitting rules for two pairs of terminal nodes
                       curr_tree,         # The current set of trees (not required if type is stump)
                       node_min_size,     # The minimum size of a node to grow
                       s,                 # probability vector to be used during the growing process
                       common_vars,       # common variables between the subsets x1 and x2
                       aux_factor_var)    # identify factor variables and
{
  
  # Call the appropriate function to get the new tree
  
  new_tree = switch(type,
                    grow = grow_tree(X, y, curr_tree, node_min_size, s, common_vars, aux_factor_var),
                    prune = prune_tree(X, y, curr_tree),
                    change = change_tree(X, y, curr_tree, node_min_size, common_vars, aux_factor_var))
  
  vars_tree = new_tree$tree_matrix[,'split_variable'] # get the split variables
  vars_tree_no_NAs = unique(vars_tree[!is.na(vars_tree)]) # remove the NAs
  # This can happen due to a previous prune step, but as soon as it happens, we set it to a stump
  if (all(vars_tree_no_NAs %in% common_vars) && length(vars_tree_no_NAs) == 1){
    new_tree = create_stump(1, y, X)[[1]]
    new_tree$ForceStump = TRUE
  }
  
  # Return the new tree
  return(new_tree)
  
} 
# Grow_tree function ------------------------------------------------------

grow_tree = function(X, y, curr_tree, node_min_size, s, common_vars, aux_factor_var) {
  
  available_values = NULL
  max_bad_trees = 10
  count_bad_trees = 0
  bad_trees = TRUE
  
  while (bad_trees ){
    
    # Set up holder for new tree
    new_tree = curr_tree
    
    # Get the list of terminal nodes
    terminal_nodes = as.numeric(which(new_tree$tree_matrix[,'terminal'] == 1))
    
    # Find terminal node sizes
    terminal_node_size = as.numeric(new_tree$tree_matrix[terminal_nodes,'node_size'])
    
    # Add two extra rows to the tree in question
    new_tree$tree_matrix = rbind(new_tree$tree_matrix,
                                 c(1, NA, NA, NA, NA, NA, NA, NA), # Make sure they're both terminal
                                 c(1, NA, NA, NA, NA, NA, NA, NA))
    
    # Choose a random terminal node to split
    node_to_split = sample(terminal_nodes, 1,
                           prob = as.integer(terminal_node_size > node_min_size)) # Choose which node to split, set prob to zero for any nodes that are too small
    
    # Choose a split variable uniformly from all columns (the first one is the intercept)
    terminal_ancestors = get_ancestors(curr_tree) # get the ancestor for all terminal nodes
    split_variable = sample(1:ncol(X), 1, prob = s)
    node_ancestors = unique(c(split_variable, terminal_ancestors[terminal_ancestors[,1] == node_to_split,2])) # covariates used in the splitting rules of the ancestor nodes + new_variable
    check_validity_new_tree = !any(unlist(lapply(aux_factor_var, function(x) all(node_ancestors %in% x)))) # check whether the new structure is valid
    
    if(check_validity_new_tree == FALSE && length(node_ancestors) == 1) {check_validity_new_tree = TRUE}
    
    # Alternatively follow BARTMachine and choose a split value using sample on the internal values of the available
    available_values = sort(unique(X[new_tree$node_indices == node_to_split,
                                     split_variable]))
    
    if(length(available_values) == 1){
      split_value = available_values[1]
    } else if (length(available_values) == 2){
      split_value = available_values[2]
    }  else {
      # split_value = sample(available_values[-c(1,length(available_values))], 1)
      split_value = resample(available_values[-c(1,length(available_values))])
    }
    
    curr_parent = new_tree$tree_matrix[node_to_split, 'parent'] # Make sure to keep the current parent in there. Will be NA if at the root node
    new_tree$tree_matrix[node_to_split,1:6] = c(0, # Now not terminal
                                                nrow(new_tree$tree_matrix) - 1, # child_left is penultimate row
                                                nrow(new_tree$tree_matrix),  # child_right is penultimate row
                                                curr_parent,
                                                split_variable,
                                                split_value)
    
    #  Fill in the parents of these two nodes
    new_tree$tree_matrix[nrow(new_tree$tree_matrix),'parent'] = node_to_split
    new_tree$tree_matrix[nrow(new_tree$tree_matrix)-1,'parent'] = node_to_split
    
    # Now call the fill function on this tree
    new_tree = fill_tree_details(new_tree, X)
    
    # Store the covariate name to use it to update the Dirichlet prior of Linero (2016).
    new_tree$var = split_variable
    
    # Check for bad tree
    if(any(as.numeric(new_tree$tree_matrix[,'node_size']) <= node_min_size) || check_validity_new_tree == FALSE) {
      count_bad_trees = count_bad_trees + 1
    } else {
      # Check whether the split variable is common to x1
      if (split_variable %in% common_vars && length(node_ancestors) == 1){ # double grow if the tree is a stump and if the split variable at the top of the tree is common x1 and x2
        s_aux = s
        s_aux[aux_factor_var[[split_variable]]] = 0 # set zero to the probability of the split variable that was just added in the tree, which is common to x1
        new_tree_double_grow = grow_tree(X, y, new_tree, node_min_size, s_aux, common_vars, aux_factor_var)
        new_tree_double_grow$var[2] = new_tree$var
        # save = new_tree$var
        new_tree_double_grow$ForceStump = FALSE
        return(new_tree_double_grow)
      }
      bad_trees = FALSE
    }
    
    if(count_bad_trees == max_bad_trees) {
      curr_tree$var = 0
      return(curr_tree)
    }
  }
  
  # Create an auxiliary variable
  new_tree$ForceStump = FALSE
  
  # Return new_tree
  return(new_tree)
  
} # End of grow_tree function

# Prune_tree function -----------------------------------------------------

prune_tree = function(X, y, curr_tree) {
  
  # Create placeholder for new tree
  new_tree = curr_tree
  
  if(nrow(new_tree$tree_matrix) == 1) { # No point in pruning a stump!
    new_tree$var = 0
    return(new_tree)
  }
  
  # Get the list of terminal nodes
  terminal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 1)
  
  # Pick a random terminal node to prune
  # ONLY PICK NODES WHERE BOTH LEFT AND RIGHT CHILD ARE TERMINAL
  bad_node_to_prune = TRUE # Assume a bad node pick
  while(bad_node_to_prune) {
    
    # Choose a random terminal node
    node_to_prune = resample(terminal_nodes)
    
    # Find the parent of this terminal node
    parent_pick = as.numeric(new_tree$tree_matrix[node_to_prune, 'parent'])
    var_pruned_nodes = as.numeric(new_tree$tree_matrix[parent_pick, 'split_variable'])
    
    # Get the two children of this parent
    child_left = as.numeric(new_tree$tree_matrix[parent_pick, 'child_left'])
    child_right = as.numeric(new_tree$tree_matrix[parent_pick, 'child_right'])
    
    # See whether either are terminal
    child_left_terminal = as.numeric(new_tree$tree_matrix[child_left, 'terminal'])
    child_right_terminal = as.numeric(new_tree$tree_matrix[child_right, 'terminal'])
    
    # If both are terminal then great
    if( (child_left_terminal == 1) & (child_right_terminal == 1) ) {
      bad_node_to_prune = FALSE # Have chosen a pair of terminal nodes so exist while loop
    }
    
  }# End of bad node to prune while loop
  
  # Delete these two rows from the tree matrix
  new_tree$tree_matrix = new_tree$tree_matrix[-c(child_left,child_right),,
                                              drop = FALSE]
  # Make this node terminal again with no children or split values
  new_tree$tree_matrix[parent_pick,c('terminal',
                                     'child_left',
                                     'child_right',
                                     'split_variable',
                                     'split_value')] = c(1, NA, NA, NA, NA)
  
  # If we're back to a stump no need to call fill_tree_details
  if(nrow(new_tree$tree_matrix) == 1) {
    new_tree$var = var_pruned_nodes
    new_tree$node_indices = rep(1, length(y))
  } else {
    # If we've removed some nodes from the middle we need to re-number all the child_left and child_right values - the parent values will still be correct
    if(node_to_prune <= nrow(new_tree$tree_matrix)) { # Only need do this if we've removed some observations from the middle of the tree matrix
      # If you're pruning any nodes which affect parent indices further down the tree then make sure to shift the parent values
      bad_parents = which(as.numeric(new_tree$tree_matrix[,'parent'])>=node_to_prune)
      # Shift them back because you have removed two rows
      new_tree$tree_matrix[bad_parents,'parent'] = as.numeric(new_tree$tree_matrix[bad_parents,'parent']) - 2
      
      for(j in node_to_prune:nrow(new_tree$tree_matrix)) {
        # Find the current parent
        curr_parent = as.numeric(new_tree$tree_matrix[j,'parent'])
        # Find both the children of this node
        curr_children = which(as.numeric(new_tree$tree_matrix[,'parent']) == curr_parent)
        # Input these children back into the parent
        new_tree$tree_matrix[curr_parent,c('child_left','child_right')] = sort(curr_children)
      } # End for loop of correcting parents and children
    } # End if statement to fill in tree details
    
    # Call the fill function on this tree
    new_tree = fill_tree_details(new_tree, X)
    
    # Store the covariate name that was used in the splitting rule of the terminal nodes that were just pruned
    new_tree$var = var_pruned_nodes
  }
  
  # Create an auxiliary variable
  new_tree$ForceStump = FALSE
  
  # Return new_tree
  return(new_tree)
  
} # End of prune_tree function

# change_tree function ----------------------------------------------------

change_tree = function(X, y, curr_tree, node_min_size, common_vars, aux_factor_var) {
  
  # Change a node means change out the split value and split variable of an internal node. Need to make sure that this does now produce a bad tree (i.e. zero terminal nodes)
  
  # If current tree is a stump nothing to change
  if(nrow(curr_tree$tree_matrix) == 1) {
    curr_tree$var = c(0, 0)
    return(curr_tree)
  }
  
  # Create a holder for the new tree
  new_tree = curr_tree
  
  # Select internal nodes which are parent of two terminals
  # internal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 0)
  terminal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 1)
  internal_parent_of_terminals = table(new_tree$tree_matrix[terminal_nodes,'parent'])
  internal_parent_of_two_terminals = which(internal_parent_of_terminals > 1)
  internal_nodes = as.numeric(names(internal_parent_of_terminals[internal_parent_of_two_terminals]))
  
  # Create a while loop to get good trees
  # Create a counter to stop after a certain number of bad trees
  max_bad_trees = 2
  count_bad_trees = 0
  bad_trees = TRUE
  while(bad_trees) {
    # Re-set the tree
    new_tree = curr_tree
    
    # choose an internal node to change
    node_to_change = resample(internal_nodes)
    
    # get the ancestors for internal nodes
    save_ancestor = get_ancestors_internal(new_tree)
    aux_internals_below_above = xtabs(split_var ~ internal + parent, save_ancestor)
    index_internals_above = which(rownames(aux_internals_below_above)==node_to_change)
    index_internals_below = which(colnames(aux_internals_below_above)==node_to_change)
    
    store_internal_above = colnames(aux_internals_below_above)[which(aux_internals_below_above[index_internals_above,] != 0)]
    store_internal_below = rownames(aux_internals_below_above)[which(aux_internals_below_above[,index_internals_below] != 0)]
    internal_nodes_above_below = as.numeric(c(store_internal_above, store_internal_below))
    split_vars_internal_nodes_above_below = unique(new_tree$tree_matrix[internal_nodes_above_below, 'split_variable'])
    
    # Get the covariate that will be changed
    var_changed_node = as.numeric(new_tree$tree_matrix[node_to_change, 'split_variable'])
    
    # Use the get_children function to get all the children of this node
    all_children = get_children(new_tree$tree_matrix, node_to_change)
    
    # Now find all the nodes which match these children
    use_node_indices = !is.na(match(new_tree$node_indices, all_children))
    
    # Create new split variable and value based on ignorance
    # then check this doesn't give a bad tree
    
    available_values = NULL
    
    new_split_variable = sample(1:ncol(X), 1)
    
    node_ancestors = unique(c(new_split_variable, split_vars_internal_nodes_above_below)) # covariates used in the splitting rules of the ancestor nodes + new_variable
    check_validity_new_tree = !any(unlist(lapply(aux_factor_var, function(x) all(node_ancestors %in% x)))) # check whether the new structure is valid
    
    if(check_validity_new_tree == FALSE && length(node_ancestors) == 1) {check_validity_new_tree = TRUE}
    
    available_values = sort(unique(X[use_node_indices,
                                     new_split_variable]))
    
    if (length(available_values) == 1){
      new_split_value = available_values[1]
      new_tree$var = c(var_changed_node, new_split_variable)
    } else if (length(available_values) == 2){
      new_split_value = available_values[2]
      new_tree$var = c(var_changed_node, new_split_variable)
    } else {
      # new_split_value = sample(available_values[-c(1,length(available_values))], 1)
      new_split_value = resample(available_values[-c(1,length(available_values))])
    }
    # Update the tree details
    new_tree$tree_matrix[node_to_change,
                         c('split_variable',
                           'split_value')] = c(new_split_variable,
                                               new_split_value)
    
    # Update the tree node indices
    new_tree = fill_tree_details(new_tree, X)
    
    # Store the covariate name that was used in the splitting rule of the terminal node that was just changed
    new_tree$var = c(new_split_variable, var_changed_node)
    
    # Check for bad tree
    node_side_aux = new_tree$tree_matrix[terminal_nodes, 'node_size']
    split_var_aux1 = new_tree$tree_matrix[1,'split_variable'] # for trees with 2 terminal nodes only
    
    if(any(as.numeric(node_side_aux) <= node_min_size) ||
       (nrow(new_tree$tree_matrix) == 3 && split_var_aux1 %in% common_vars) || check_validity_new_tree == FALSE) {
      count_bad_trees = count_bad_trees + 1
    } else {
      bad_trees = FALSE
    }
    if(count_bad_trees == max_bad_trees){
      curr_tree$var = c(0, 0)
      return(curr_tree)
    }
    
  } # end of while loop
  
  # Create an auxiliary variable
  new_tree$ForceStump = FALSE
  
  # Return new_tree
  return(new_tree)
  
} # End of change_tree function

# swap_tree function ------------------------------------------------------

swap_tree = function(X, y, curr_tree, node_min_size) {
  
  # Swap takes two neighbouring internal nodes and swaps around their split values and variables
  
  # If current tree is a stump nothing to change
  if(nrow(curr_tree$tree_matrix) == 1) return(curr_tree)
  
  # Create a holder for the new tree
  new_tree = curr_tree
  
  # Need to get the internal nodes
  internal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 0)
  terminal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 1)
  
  # If less than 3 internal nodes return curr_tree
  if(length(internal_nodes) < 3) return(curr_tree)
  
  # Find pairs of neighbouring internal nodes
  parent_of_internal = as.numeric(new_tree$tree_matrix[internal_nodes,'parent'])
  pairs_of_internal = cbind(internal_nodes, parent_of_internal)[-1,]
  
  # Create a while loop to get good trees
  # Create a counter to stop after a certain number of bad trees
  max_bad_trees = 2
  count_bad_trees = 0
  bad_trees = TRUE
  while(bad_trees) {
    # Re-set the tree
    new_tree = curr_tree
    
    # Pick a random pair
    nodes_to_swap = sample(1:nrow(pairs_of_internal), 1)
    
    # Get the split variables and values for this pair
    swap_1_parts = as.numeric(new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,1],
                                                   c('split_variable', 'split_value')])
    swap_2_parts = as.numeric(new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,2],
                                                   c('split_variable', 'split_value')])
    
    # Update the tree details - swap them over
    new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,1],
                         c('split_variable',
                           'split_value')] = swap_2_parts
    new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,2],
                         c('split_variable',
                           'split_value')] = swap_1_parts
    
    # Update the tree node indices
    new_tree = fill_tree_details(new_tree, X)
    
    # Check for bad tree
    if(any(as.numeric(new_tree$tree_matrix[terminal_nodes, 'node_size']) <= node_min_size)) {
      count_bad_trees = count_bad_trees + 1
    } else {
      bad_trees = FALSE
    }
    if(count_bad_trees == max_bad_trees) return(curr_tree)
    
  } # end of while loop
  
  # Return new_tree
  return(new_tree)
  
} # End of swap_tree function


# Get predictions给定树和x，给定预测值 ---------------------------------------------------------
get_predictions = function(trees, X, single_tree = FALSE) {
  
  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]
  
  # Normally trees will be a list of lists but just in case
  if(single_tree) {
    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(trees, X)$node_indices
      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        predictions[curr_X_node_indices == unique_node_indices[i]] =
          trees$tree_matrix[unique_node_indices[i], 'mu']
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {
    # Do a recursive call to the function
    partial_trees = trees
    partial_trees[[1]] = NULL # Blank out that element of the list
    predictions = get_predictions(trees[[1]], X, single_tree = TRUE)  +
      get_predictions(partial_trees, X,
                      single_tree = length(partial_trees) == 1)
    #single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }
  
  return(predictions)
}


# Get tree priors ---------------------------------------------------------

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

sample_move = function(curr_tree, i, nburn){
  
  if (nrow(curr_tree$tree_matrix) == 1 || i < max(floor(0.1*nburn), 10)) {
    type = 'grow'
  } else {
    type = sample(c('grow', 'prune', 'change'), 1)
  }
  return(type)
}

get_ancestors = function(tree){
  
  save_ancestor = NULL
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  
  if(nrow(tree$tree_matrix) == 1) {
    save_ancestor = cbind(terminal = NULL,
                          ancestor = NULL)
  } else {
    for (k in seq_len(length(which_terminal))) {
      get_parent = as.numeric(as.character(tree$tree_matrix[which_terminal[k], 'parent'])) # get the 1st parent
      get_split_var = as.character(tree$tree_matrix[get_parent, 'split_variable']) # then, get the covariate associated to the row of the parent
      
      save_ancestor = rbind(save_ancestor,
                            cbind(terminal = which_terminal[k],
                                  # parent   = get_parent,
                                  ancestor = get_split_var))
      while (get_parent > 1){
        get_parent = as.numeric(as.character(tree$tree_matrix[get_parent,'parent'])) # then, get the subsequent parent
        get_split_var = as.character(tree$tree_matrix[get_parent, 'split_variable']) # then, get the covariate associated to the row of the new parent
        save_ancestor = rbind(save_ancestor,
                              cbind(terminal = which_terminal[k],
                                    # parent   = get_parent,
                                    ancestor = get_split_var))
      }
    }
    save_ancestor = unique(save_ancestor) # remove duplicates
    save_ancestor = save_ancestor[order(save_ancestor[,1], save_ancestor[,2]),] # sort by terminal and ancestor
  }
  
  return(save_ancestor)
}

resample <- function(x, ...) x[sample.int(length(x), size=1), ...]

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

get_ancestors_internal = function(tree){
  save_ancestor = NULL
  tree = tree$tree_matrix
  which_internal = which(tree[,'terminal'] == 0)
  
  
  if(nrow(tree) == 1) {
    save_ancestor = cbind(internal = NULL,
                          ancestor = NULL)
  } else {
    for (k in length(which_internal):1) {
      internal_node = which_internal[k]
      parent = tree[internal_node, 'parent']
      get_split_var = tree[internal_node, 'split_variable']
      
      save_ancestor = rbind(save_ancestor,
                            cbind(internal = internal_node,
                                  parent   = parent,
                                  split_var = get_split_var))
      while (is.na(parent) == FALSE && parent > 0) {
        get_split_var = tree[parent, 'split_variable']
        parent = tree[parent, 'parent']
        save_ancestor = rbind(save_ancestor,
                              cbind(internal = internal_node,
                                    parent   = parent,
                                    split_var = get_split_var))
      }
    }
  }
  return(save_ancestor[,,drop=FALSE])
  
}

get_children = function(tree_mat, parent) {
  # Create a holder for the children
  all_children = NULL
  if(as.numeric(tree_mat[parent,'terminal']) == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left = as.numeric(tree_mat[parent, 'child_left'])
    curr_child_right = as.numeric(tree_mat[parent, 'child_right'])
    # Return the children and also the children of the children recursively
    return(c(all_children,
             get_children(tree_mat,curr_child_left),
             get_children(tree_mat,curr_child_right)))
  }
}