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
## Posterior sampling for individual state (Ohio/New York)
ohio_beta <- function(x, y, m, inits, d, n.chains=2, n.adapt=500, n.iter=10000){
  dat <- list("x"=x, "y"=y, "m" = m, "d" = d)
  cat("model
      {
        for(i in 1:m){
          p[i] <- exp(inprod(beta,x[i,]))/(1+exp(inprod(beta,x[i,])))
          y[i] ~ dbern(p[i])
        }
        
        for(k in 1:d){
          for(j in 1:d){
            sigma[k,j] <- equals(k,j)
          }
        }
  
        mu <- rep(0,d)
        beta ~ dmnorm(mu,sigma)
      }",file="ohio_jags.txt")
  jags.m <- jags.model(file="ohio_jags.txt", data=dat, inits = inits,
                       n.chains=n.chains, n.adapt=n.adapt)
  params <- "beta"
  samps <- coda.samples(jags.m, params, n.iter=n.iter) 
  return(samps)
}
## Log likelihood function
log_L_func <-function(beta, X_ny, y_ny, m){
  res = 0
  for(i in 1:m){
    res = res + y_ny[i]*(X_ny[i,]%*%t(beta))-log(1+exp(X_ny[i,]%*%t(beta)))
  }
  return(res)
}
## Graph-enabled MCMC
### Parameters: n.chain: number of chains
###             chain.length: length of one chain
###             k, h, rho: algorithm parameters
###             Theta: prior samples 
###             X_ny, y_ny: observation matrix 
mcmc.sim <- function(n.chain=3, chain.length, k, h, rho, Theta, X_ny, y_ny){
  
  B = dim(Theta)[1]
  d = dim(Theta)[2]
  m = dim(X_ny)[1]
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
  idx <- !duplicated(t(apply(graph_edgelist,  1, sort))) 
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
    A[i, nbh.list[[i]]] <- 1
  }
  DN = NULL
  for(i in 1:B){
    DN = c(DN, length(nbh.list[[i]]))
  }
  
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
      temp <- c(log_L_func(beta_tilde, X_ny, y_ny, m), log_L_func(alpha_tilde, X_ny, y_ny, m),
                rho/B+indicator*(1-rho)/DN[beta_idx], rho/B+indicator*(1-rho)/DN[alpha_idx])
      
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
    res_list <- list("init_time"=init.time, "run_time" = run.time,
                     "alpha_df" = alpha_tilde_vec)
    return_list[[c]] <- res_list
    print(paste0("Graph-enabled MCMC: Chain ",c))
  }
  return(return_list)
}
## Metropolis random walk
### KDE evaluation
pi_hat_func <- function(theta, X, h){
  res <- 0
  B <- nrow(X)
  for(i in 1:B){
    res <- res + dnorm(norm(as.matrix(theta-X[i,]), type="f")/h, 0, 1)
  }
  return(res/(B*(h^d)))
}
### Parameters: n.chain: number of chains
###             chain.length: length of one chain
###             sigma_eps, h: algorithm parameters
###             Theta: prior samples 
###             X_ny, y_ny: observation matrix 
MRW <- function(n.chain = 3, chain.length, sigma_eps, h, Theta, X_ny, y_ny){
  
  B = dim(Theta)[1]
  d = dim(Theta)[2]
  m = dim(X_ny)[1]
  return_list = NULL
  
  run.time <- rep(0, chain.length)
  
  for(c in 1:n.chain){
    
    start.time <- Sys.time()
    
    alpha <- t(as.matrix(Theta[sample(1:B,1),]))
    alpha_list <- alpha
    
    run.time[1] <- Sys.time() - start.time
    t <- 1
    
    while(t<chain.length){
      start.time <- Sys.time()
      eps <- rmvnorm(1, mean = rep(0,d), sigma = diag(sigma_eps^2,d))
      alpha_new <- alpha + eps
      temp <- c(pi_hat_func(alpha_new,Theta,h), log_L_func(alpha_new, X_ny, y_ny, m),
                pi_hat_func(alpha,Theta,h), log_L_func(alpha, X_ny, y_ny, m))
      log_MH_ratio <- log(temp[1])+temp[2]-log(temp[3])-temp[4]
      MH_ratio <- exp(log_MH_ratio)
      if(runif(1,0,1) < MH_ratio){
        alpha <- alpha_new
      }
      alpha_list = rbind(alpha_list, alpha)
      t <- t+1
      run.time[t] <- Sys.time() - start.time
    }
    res_list <- list("run_time" = run.time,
                     "alpha_df" = alpha_list)
    return_list[[c]] <- res_list
    print(paste0("Metropolis random walk: Chain ",c))
  }
  return(return_list)
}

# Generating the data
Generate_data <- function(m, d, mu1, mu2){
  mu_0 <- rep(0, d)
  sigma_0 <- diag(1, d, d)
  beta_0 <- rmvnorm(1, mean = mu_0, sigma = sigma_0)
  beta <- as.matrix(beta_0)
  #Ohio data generation
  mu_1 <- rep(mu1, d)
  sigma_1 <- diag(1, d, d)
  X_ohio <- rmvnorm(m, mean = mu_1, sigma = sigma_1)
  log_term <- exp(X_ohio%*%t(beta))
  y_ohio <- rep(NA, m)
  for(i in 1:m){
    y_ohio[i] <- rbinom(1, 1, prob = log_term[i]/(1+log_term[i]))
  }
  #New York data generation
  mu_2 <- rep(mu2, d)
  sigma_2 <- diag(1, d, d)
  X_ny <- rmvnorm(m, mean = mu_2, sigma = sigma_2)
  log_term <- exp(X_ny%*%t(beta))
  y_ny <- rep(NA, m)
  for(i in 1:m){
    y_ny[i] <- rbinom(1, 1, prob = log_term[i]/(1+log_term[i]))
  }
  dic_name <- paste0("./OUD/d=", d, ",")
  saveRDS(X_ohio, paste0(dic_name,"X_ohio.rds"))
  saveRDS(y_ohio, paste0(dic_name,"y_ohio.rds"))
  saveRDS(X_ny, paste0(dic_name,"X_ny.rds"))
  saveRDS(y_ny, paste0(dic_name,"y_ny.rds"))
  saveRDS(beta_0, paste0(dic_name,"beta0.rds"))
}

# Posterior sampling
Posterior_sampling <- function(m, d, n.iter = 40000){
  dic_name <- paste0("./OUD/d=", d, ",")
  X_ohio <- readRDS(paste0(dic_name,"X_ohio.rds"))
  y_ohio <- readRDS(paste0(dic_name,"y_ohio.rds"))
  X_ny <- readRDS(paste0(dic_name,"X_ny.rds"))
  y_ny <- readRDS(paste0(dic_name,"y_ny.rds"))
  #Ohio posterior
  samps_ohio <- ohio_beta(x=X_ohio, y=y_ohio, m=m, d=d, inits=NULL,
                          n.chains=3, n.adapt=5000, n.iter=n.iter)
  print(gelman.diag(samps_ohio))
  #New York posterior with uninformative prior
  samps_ny <- ohio_beta(x=X_ny, y=y_ny, m=m, d=d, inits=NULL,
                        n.chains=3, n.adapt=5000, n.iter=n.iter)
  print(gelman.diag(samps_ny))
  #New York posterior with informative prior (true)
  samps_true <- ohio_beta(x=rbind(X_ohio,X_ny), y=c(y_ohio,y_ny), m=2*m, d=d, inits=NULL, 
                          n.chains=3, n.adapt=5000, n.iter=n.iter)
  print(gelman.diag(samps_true))
  
  saveRDS(samps_ohio,file=paste0(dic_name,"samps_ohio.rds"))
  saveRDS(samps_ny,file=paste0(dic_name,"samps_ny.rds"))
  saveRDS(samps_true,file=paste0(dic_name,"samps_true.rds"))
}

# Simulations using the graph-enabled MCMC algorithm and the Metropolis random walk
Simulations <- function(B, d, n.chain, chain.length, k, h, rho, sigma_eps){
  dic_name <- paste0("./OUD/d=", d, ",")
  X_ny <- readRDS(paste0(dic_name,"X_ny.rds"))
  y_ny <- readRDS(paste0(dic_name,"y_ny.rds"))
  samps_ohio <- readRDS(paste0(dic_name,"samps_ohio.rds"))
  n.iter = dim(samps_ohio[[1]])[1]
  Theta_ohio <- samps_ohio[[1]][(n.iter-B+1):n.iter,]
  dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
  # Graph-enabled MCMC
  result <- mcmc.sim(n.chain=n.chain, chain.length=chain.length, k=k, h=h, rho=rho, Theta=Theta_ohio, X_ny=X_ny, y_ny=y_ny)
  list.save(result, paste0(dic_name_B,"result.rdata"))
  # Metropolis random walk
  res_mrw <- MRW(n.chain=n.chain, chain.length=chain.length, sigma_eps=sigma_eps, h=h, Theta=Theta_ohio, X_ny=X_ny, y_ny=y_ny)
  list.save(res_mrw, paste0(dic_name_B,"result_mrw.rdata"))
}

# Simulations
## Parameters
n.chain = 3
chain.length = 10000
m = 1500 # Number of data points in Ohio/New York
mu1 = 1
mu2 = -1
h = 0.04
rho = 0.5
sigma_eps = 0.02 
n.iter = 40000

## Data generation and posterior sampling
set.seed(0)
d = 2
Generate_data(m=m, d=d, mu1=mu1, mu2=mu2)
Posterior_sampling(m=m, d=d, n.iter=n.iter)
d = 6
Generate_data(m=m, d=d, mu1=mu1, mu2=mu2)
Posterior_sampling(m=m, d=d, n.iter=n.iter)
d = 10
Generate_data(m=m, d=d, mu1=mu1, mu2=mu2)
Posterior_sampling(m=m, d=d, n.iter=n.iter)

## Graph-enabled MCMC and Metropolis random walk
### B=10000; d=2, 6, 10
B = 10000 # Number of prior draws
k = ceiling(sqrt(B))

d = 2
set.seed(0)
Simulations(B=B, d=d, n.chain=n.chain, chain.length=chain.length, k=k, h=h, rho=rho, sigma_eps=sigma_eps)

d = 6
set.seed(0)
Simulations(B=B, d=d, n.chain=n.chain, chain.length=chain.length, k=k, h=h, rho=rho, sigma_eps=sigma_eps)

d = 10
set.seed(0)
Simulations(B=B, d=d, n.chain=n.chain, chain.length=chain.length, k=k, h=h, rho=rho, sigma_eps=sigma_eps)

### d=6; B=1000, 2500, 5000, 15000, 20000
d = 6

B <- 1000
k <- ceiling(sqrt(B))
set.seed(0)
Simulations(B=B, d=d, n.chain=n.chain, chain.length=chain.length, k=k, h=h, rho=rho, sigma_eps=sigma_eps)

B <- 2500 
k <- ceiling(sqrt(B))
set.seed(0)
Simulations(B=B, d=d, n.chain=n.chain, chain.length=chain.length, k=k, h=h, rho=rho, sigma_eps=sigma_eps)

B <- 5000 
k <- ceiling(sqrt(B))
set.seed(0)
Simulations(B=B, d=d, n.chain=n.chain, chain.length=chain.length, k=k, h=h, rho=rho, sigma_eps=sigma_eps)

B <- 15000 
k <- ceiling(sqrt(B))
set.seed(0)
Simulations(B=B, d=d, n.chain=n.chain, chain.length=chain.length, k=k, h=h, rho=rho, sigma_eps=sigma_eps)

B <- 20000 
k <- ceiling(sqrt(B))
set.seed(0)
Simulations(B=B, d=d, n.chain=n.chain, chain.length=chain.length, k=k, h=h, rho=rho, sigma_eps=sigma_eps)

