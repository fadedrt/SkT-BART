# SkT-BART
This repository contains code to implement the methodology described in the paper “Bayesian Adve Regression Trees for Responses with Heavy Tails and Skewness”. 
The tree structure and framework in this code are adapted from the primary functions of [CSP-BART](https://github.com/ebprado/CSP-BART/tree/main).

Model comparisons are performed using the [BART](https://github.com/ebprado/CSP-BART/blob/main/cspbart/R/bart.R) module from CSP-BART and the implementation of [skewBART](https://github.com/Seungha-Um/skewBART?tab=readme-ov-file).

## Example
```r
library(skewt)
library(zeallot)
library(rust)
sim_fried <- function(N, P, v, sigma,gamma) {
  #lambda <- alpha * sigma/sqrt(1+alpha^2)
  #tau <- sigma/sqrt(1+alpha^2)
  X <- matrix(runif(N * P), nrow = N)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
  #X[,1] * X[,2] + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
  #Z <- abs(rnorm(N, mean=0, sd=1) )
  Y <- mu +sigma*rskt(N, v, gamma)
  #EY <- mu + lambda * sqrt(2/pi)
  return(list(X = X, Y = Y, mu = mu, v=v,gamma))
}
prediction <- function(x,y,test_x,test_y,npost,trees){
  y_mean = mean(y)
  y_sd = sd(y)
  prediction_train= matrix(NA,ncol = npost, nrow = nrow(x))
  prediction_test=matrix(NA, ncol = npost,nrow=nrow(test_x))
  for (i in 1:npost){
    prediction_train[,i] <- get_predictions(trees[[i]],x,single_tree = FALSE)
  }
  for (i in 1:npost){
    prediction_test[,i] <- get_predictions(trees[[i]],test_x,single_tree = FALSE)
  }
  
  f_hat_train      <- prediction_train * y_sd + y_mean
  f_hat_test       <- prediction_test * y_sd + y_mean
  f_hat_train_mean <- rowMeans(f_hat_train)
  f_hat_test_mean  <- rowMeans(f_hat_test)
  return(list(f_hat_train = f_hat_train, f_hat_test = f_hat_test, 
              f_hat_train_mean = f_hat_train_mean, 
              f_hat_test_mean = f_hat_test_mean))
}

# Calculation of LPML values 
calculate_lpml <- function(y, x, npost=100, trees, sigma2_store, gamma2_store, v_store) {
  y_mean <- mean(y)
  y_sd <- sd(y)
  n <- length(y)
  
  # 每个 MCMC 抽样下的预测值矩阵 (n x npost)
  prediction_trees <- matrix(NA, nrow = n, ncol = npost)
  for (j in 1:npost) {
    prediction_trees[, j] <- get_predictions(trees[[j]], x, single_tree = FALSE) * y_sd + y_mean
  }
  
  # 初始化 CPO 存储
  cpo <- numeric(n)
  for (i in 1:n) {
    predictive_densities <- numeric(npost)
    for (j in 1:npost) {
      mu_ij <- prediction_trees[i, j]
      sigma_j <- sqrt(sigma2_store[j])
      gamma_j <- sqrt(gamma2_store[j])
      v_j <- v_store[j]
      
      z <- (y[i] - mu_ij) / sigma_j
      predictive_densities[j] <- dskt(z, df = v_j, gamma = gamma_j) / sigma_j
    }
    
    predictive_densities <- pmax(predictive_densities, 1e-8)  # 避免数值下溢
    cpo[i] <- 1 / mean(1 / predictive_densities)  # 调和平均
  }
  
  # 计算 LPML
  lpml <- sum(log(cpo))
  return(lpml)
}

set.seed(123)
c(x,y,mu,v,gamma) %<-% sim_fried(250, 5, 1,1,1)
set.seed(123)
c(test_x,test_y,test_mu,v,gamma)  %<-% sim_fried(100, 5,1,1,1)
sktbart_fit = sktbart(y,x,ntrees =50,nburn=2500,npost =2500,lambda1 = runif(length(y),min=0,max = 1)) 
v_hat = mean(sktbart_fit$v)
sigma_hat = sqrt(mean(sktbart_fit$sigma2))
gamma_hat = sqrt(mean(sktbart_fit$gamma2))

if (v_hat > 1) {
  r = sigma_hat * rskt(10000, v_hat, gamma_hat)
} else if (v_hat > 0 && v_hat <= 1) {
  # When the degrees of freedom are less than or equal to 1, 
  # the expectation of the t-distribution does not exist. 
  # Therefore, for convenience, we set the adjustment factor 
  # of the error term to 0 when v_hat <= 1.
  r = 0
}

c(f_hat_train, f_hat_test, f_hat_train_mean, f_hat_test_mean) %<-% 
  prediction(x, y, test_x, test_y, npost = 2500, trees = sktbart_fit)

cat("sktBART的RMSE为:", round(sqrt(sum((f_hat_test_mean - test_y)^2)/length(test_y)), 2), "\n")
cat("sktBART的MAE为:", round(sum(abs(f_hat_test_mean - test_y))/length(test_y), 2), "\n")
cat("sktBART的MAPE为:", round(mean(abs(f_hat_test_mean - test_y) / abs(test_y)) * 100, 2), "\n")
cat("sktBART的LPML为:"calculate_lpml(y=y,x=x,npost=2500,trees=sktbart_fit$trees,sigma2_store = sktbart_fit$sigma2,gamma2_store=sktbart_fit$gamma2,v_store = sktbart_fit$v))
