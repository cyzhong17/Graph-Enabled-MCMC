# Loading packages
library(truncnorm)
library(igraph)
library(mvtnorm)
library(animation)
library(ndtv)
library(expm)
library(MVN)
library(coda)
library(mcmcse)
library(rjags)
library(gtools)
library(MASS)
library(survival)
library(gTests)
library(ggplot2)
library(e1071)
library(energy)
library(kernlab)
library(ICSNP)
library(rlist)
library(kableExtra)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(transport)
library(reshape2)

# Functions
## Log likelihood for Gaussian observations
### beta: parameter vector
log_L_func_G <-function(beta, tau, n, X){
  res = 0
  for(i in 1:n){
    res = res - sum((X[i,]-beta)^2)
  }
  res = res/(2*tau^2) - n*d*(log(2*pi)/2+log(tau))
  return(res)
}
## Graph-enabled MCMC
### Parameters: n.chain: number of chains
###             chain.length: length of one chain
###             k, h, rho: algorithm parameters
###             Theta: prior samples 
###             X: observation matrix 
mcmc.sim_G <- function(n.chain=3, chain.length, k, h, rho, Theta, X, tau){
  
  B = dim(Theta)[1]
  d = dim(Theta)[2]
  n = dim(X)[1]
  return_list = NULL
  
  start.init = Sys.time()
  
  dist_vec <- dist(Theta, method = "euclidean", diag = T)
  dist_mat <- as.matrix(dist_vec, B, B)
  diag(dist_mat) <- Inf
  graph_edgelist <- cbind(rep(1,k), order(dist_mat[1,])[1:k])
  
  for(i in 2:B){
    graph_edgelist <- rbind(graph_edgelist,
                            cbind(rep(i,k), order(dist_mat[i,])[1:k]))
  }
  
  idx <- !duplicated(t(apply(graph_edgelist, 1, sort))) 
  
  graph_edgelist <- t(apply(graph_edgelist, 1, as.character))
  
  network <- graph_from_data_frame(d=graph_edgelist[idx,], directed=F) 
  
  nbh.list <- neighbors(network, v=as.character(1))
  nbh.list <- as.numeric(names(nbh.list))
  nbh.list <- list(nbh.list)
  for(i in 2:B){
    nbh.i <- neighbors(network, v=as.character(i))
    nbh.i <- as.numeric(names(nbh.i))
    nbh.i <- list(nbh.i)
    nbh.list <- append(nbh.list, nbh.i)
  }
  
  A <- matrix(0, nrow=B, ncol=B)
  for(i in 1:B){
    A[i,nbh.list[[i]]] <- 1
  }
  DN <- rowSums(A)
  
  init.time = Sys.time() - start.init
  
  run.time = rep(0, chain.length)
  
  for(c in 1:n.chain){
    
    start.time <- Sys.time()
    
    alpha_idx <- sample(1:B, 1)
    alpha_tilde <- rmvnorm(1, mean = Theta[alpha_idx,], sigma = diag(h^2,d))
    alpha_tilde_vec <- alpha_tilde
    
    run.time[1] = Sys.time() - start.time
    
    t <- 1
    
    while(t<chain.length){
      
      start.time <- Sys.time()
      
      ind <- rbinom(1, 1, rho)
      if(ind == 0){
        beta_idx <- sample(nbh.list[[alpha_idx]], 1)
      } else {
        beta_idx <- sample(1:B, 1)
      }
      eps <- rmvnorm(1, mean = rep(0,d), sigma=diag(h^2,d))
      beta_tilde <- Theta[beta_idx,] + eps
      
      indicator <- A[alpha_idx, beta_idx]
      temp <- c(log_L_func_G(beta_tilde,tau,n,X), log_L_func_G(alpha_tilde,tau,n,X),
                rho/B+indicator*(1-rho)/DN[beta_idx],
                rho/B+indicator*(1-rho)/DN[alpha_idx])
      log_MH_ratio <- log(temp[3])+temp[1]-log(temp[4])-temp[2]
      MH_ratio <- exp(log_MH_ratio)
      
      if(runif(1,0,1) < MH_ratio){
        alpha_idx = beta_idx
        alpha_tilde = beta_tilde
        alpha_tilde_vec <- rbind(alpha_tilde_vec, beta_tilde)
      }else{
        alpha_tilde_vec <- rbind(alpha_tilde_vec, alpha_tilde)
      }
      t <- t+1
      
      run.time[t] = Sys.time() - start.time
    }
    res_list <- list("init_time" = init.time, "run_time" = run.time,
                     "alpha_df" = alpha_tilde_vec)
    return_list[[c]] <- res_list
    print(paste0("Graph-enabled MCMC: Chain ",c))
  }
  return(return_list)
}
## Metropolis random walk
pi_hat_func_G <- function(theta, X, h){
  res <- 0
  B <- nrow(X)
  for(i in 1:B){
    res <- res + dnorm(norm(as.matrix(theta-X[i,]), type="f")/h,0,1)
  }
  return(res/(B*(h^d)))
}
### Parameters: n.chain: number of chains
###             chain.length: length of one chain
###             sigma_eps, h: algorithm parameters
###             Theta: prior samples 
###             X: observation matrix 
MRW_G <- function(n.chain = 3, chain.length = 10000, sigma_eps, h, Theta, X, tau){
  
  B = dim(Theta)[1]
  d = dim(Theta)[2]
  n = dim(X)[1]
  return_list = NULL
  
  run.time <- rep(0, chain.length)
  
  for(c in 1:n.chain){
    
    start.time <- Sys.time()
    
    alpha <- Theta[sample(1:B,1),]
    alpha_list <- alpha
    run.time[1] <- Sys.time() - start.time
    
    t <- 1
    
    while(t<chain.length){
      start.time <- Sys.time()
      eps <- rmvnorm(1, mean = rep(0,d), sigma = diag(sigma_eps^2,d))
      alpha_new <- alpha + eps
      temp <- c(pi_hat_func_G(alpha_new,Theta,h),log_L_func_G(alpha_new,tau,n,X),
                pi_hat_func_G(alpha,Theta,h),log_L_func_G(alpha,tau,n,X))
      
      log_MH_ratio <- log(temp[1])+temp[2]-log(temp[3])-temp[4]
      MH_ratio <- exp(log_MH_ratio)
      if(runif(1,0,1) < MH_ratio){
        alpha <- alpha_new
      }
      alpha_list = rbind(alpha_list, alpha)
      t <- t+1
      run.time[t] <- Sys.time() - start.time
    }
    res_list <- list("run_time" = run.time, "alpha_df" = alpha_list)
    return_list[[c]] <- res_list
    print(paste0("Metropolis random walk: Chain ",c))
  }
  return(return_list)
}

# Simulations
## Setting the parameters
sigma = 1
tau = 2
n = 10
h = 1
mu1 = c(4,0)
mu2 = c(-4,0)
mu3 = c(0,4)
B = 100

## Generating the data
set.seed(0)
I = rmultinom(1, 1, c(1/3,1/3,1/3))
if(I[1] == 1){
  theta = mu1 + rnorm(2, mean=0, sd=sigma)
} else if (I[2] == 1){
  theta = mu2 + rnorm(2, mean=0, sd=sigma)
} else {
  theta = mu3 + rnorm(2, mean=0, sd=sigma)
}
X = matrix(0, nrow=n, ncol=2)
for(i in 1:n){
  X[i,] = theta + rnorm(2, mean=0, sd=tau)
}
Theta0 = matrix(0, nrow=B, ncol=2)
for(i in 1:B){
  I = rmultinom(1, 1, c(1/3,1/3,1/3))
  if(I[1] == 1){
    Theta0[i,] = mu1 + rnorm(2, mean=0, sd=sigma)
  } else if (I[2] == 1){
    Theta0[i,] = mu2 + rnorm(2, mean=0, sd=sigma)
  } else {
    Theta0[i,] = mu3 + rnorm(2, mean=0, sd=sigma)
  }
}
## Save the data
dic_name <- "/Users/chenyangzhong/Desktop/Data_paper/GM/"
saveRDS(X, paste0(dic_name,"X.rds"))
saveRDS(Theta0, paste0(dic_name,"Theta0.rds"))
saveRDS(theta, paste0(dic_name,"theta.rds"))

## Generating the MCMC samples
### Read data
dic_name <- "/Users/chenyangzhong/Desktop/Data_paper/GM/"
X <- readRDS(paste0(dic_name,"X.rds"))
Theta0 <- readRDS(paste0(dic_name,"Theta0.rds"))
theta <- readRDS(paste0(dic_name,"theta.rds"))
### Run the chains
n.chain <- 3
chain.length <- 10000
k <- ceiling(sqrt(B))
rho <- 0.5
sigma_eps <- 0.5
result <- mcmc.sim_G(n.chain=n.chain, chain.length = chain.length, k=k, h=h, rho=rho, Theta=Theta0, X=X, tau=tau)
list.save(result, paste0(dic_name,"result.rdata"))
res_mrw <- MRW_G(n.chain=n.chain, chain.length=chain.length, sigma_eps=sigma_eps, h=h, Theta=Theta0, X=X, tau=tau)
list.save(res_mrw, paste0(dic_name,"result_mrw.rdata"))

## True posterior distribution
pi_func_true <- function(theta){
  res <- 1/3*dmvnorm(theta, mean = mu1, sigma = sigma^2*diag(1,2))
  + 1/3*dmvnorm(theta, mean = mu2, sigma = sigma^2*diag(1,2))
  + 1/3*dmvnorm(theta, mean = mu3, sigma = sigma^2*diag(1,2))
  return(res)
}
MRW_true <- function(n.chain = 3, chain.length = 10000, sigma_eps, X, tau){
  
  return_list = NULL
  d = dim(X)[2]
  n = dim(X)[1]
  
  for(c in 1:n.chain){
    
    alpha <- colMeans(X)
    alpha_list <- alpha
    
    t <- 1
    while(t<chain.length){
      eps <- rmvnorm(1, mean = rep(0,d), sigma = diag(sigma_eps^2,d))
      alpha_new <- alpha + eps
      temp <- c(pi_func_true(alpha_new),log_L_func_G(alpha_new,tau,n,X),
                pi_func_true(alpha),log_L_func_G(alpha,tau,n,X))
      log_MH_ratio <- log(temp[1])+temp[2]-log(temp[3])-temp[4]
      MH_ratio <- exp(log_MH_ratio)
      if(runif(1,0,1) < MH_ratio){
        alpha <- alpha_new
      }
      alpha_list = rbind(alpha_list, alpha)
      t <- t+1
      print(t)
    }
    res_list <- list("alpha_df" = alpha_list)
    return_list[[c]] <- res_list
    print(c)
  }
  return(return_list)
}
result_true <- MRW_true(n.chain=n.chain, chain.length=chain.length, sigma_eps=sigma_eps, X=X, tau=tau)
list.save(result_true, paste0(dic_name,"result_true.rdata"))


