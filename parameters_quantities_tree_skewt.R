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