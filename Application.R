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
library(gridExtra)

set.seed(0)

# Downloading data and generating posterior draws for the reference and current studies
## copy links for raw dataset files
gallup_coop <- url("https://raw.githubusercontent.com/gal-zz-lup/NGS2/master/cooperation_exp1_FINAL.csv")
gallup_rewire <- url("https://raw.githubusercontent.com/gal-zz-lup/NGS2/master/rewire_exp1.csv")
temp <- tempfile()
download.file("http://davidrand-cooperation.com/s/Rand2011PNAS_data_and_code-pi6b.zip",
              temp, mode="wb")
gallup1 <- read.table(unz(temp, "Rand2011PNAS_cooperation_data.txt"),
                      sep="\t", skip=0, header=T)
unlink(temp)
gallup2 <- read.csv(gallup_coop, header = TRUE, sep = ',')
names(gallup2)[names(gallup2) %in%
                 c("round","pid","action")] <- c("round_num",
                                                 "playerid","decision..0.D.1.C.")
gallup2$sessionF <- as.factor(gallup2$session)
gallup2$sessionnum <- unlist(lapply(1:nrow(gallup2), function(i){
  which(gallup2$sessionF[i]==levels(gallup2$sessionF)) }))

public_game <- function(dataset, inits, n.chains=2, n.adapt=500, n.iter=10000){
  dat <- list("x"=dataset[,c("fluid_dummy","round_num")], "y"=dataset[,"decision..0.D.1.C."],
              'm'=nrow(dataset))
  cat("model
      {
        for(i in 1:m){
          p[i] <- 1/( 1 + exp(-(beta[1]+beta[2]*x[i,1]+beta[3]*x[i,2]+beta[4]*x[i,1]*x[i,2])) )
          y[i] ~ dbern(p[i])
        }
        
        for(k in 1:4){
          for(j in 1:4){
            pre[k,j] <- equals(k,j)*precision[k]
          }
        }
        
        precision <- c(1/(2.5^2),1/(5.3^2),1/(0.64),1)
        mu <- rep(0,4)
        beta ~ dmnorm(mu,pre)
      }",file = "public_game.txt")
  jags.m <- jags.model(file="public_game.txt", data=dat, inits=inits,
                       n.chains=n.chains, n.adapt=n.adapt)
  params <- "beta"
  samps <- coda.samples(jags.m, params, n.iter=n.iter) 
  return(samps)
}
n.iter <- 20000
res_pre_pilot <- public_game(dataset=gallup1, inits = NULL,
                             n.chains=3, n.adapt=10000, n.iter=n.iter)
res_pre_exp <- public_game(dataset=gallup2, inits = NULL,
                           n.chains=3, n.adapt=10000, n.iter=n.iter)
gelman.diag(res_pre_pilot)
summary(res_pre_pilot)
gelman.diag(res_pre_exp)
summary(res_pre_exp)

## Saving the results
dic_name <- "/Users/chenyangzhong/Desktop/Data_paper/Application/"
saveRDS(res_pre_pilot, paste0(dic_name,"Posterior_pilot.rds"))
saveRDS(res_pre_exp, paste0(dic_name,"Posterior_exp.rds"))

# Visualizing posterior draws from the reference and current studies
d <- 4
B <- 5000 #number of posterior draws from the reference study
n.iter <- 20000
dic_name <- "/Users/chenyangzhong/Desktop/Data_paper/Application/"
res_pre_pilot <- readRDS(paste0(dic_name,"Posterior_pilot.rds"))
res_pre_exp <- readRDS(paste0(dic_name,"Posterior_exp.rds"))
Theta_pilot <- res_pre_pilot[[1]][(n.iter-B+1):n.iter,]
Theta_exp <- res_pre_exp[[1]][(n.iter-B+1):n.iter,]

df_pilot <- data.frame(Theta_pilot)
df_exp <- data.frame(Theta_exp)
df_pilot_long <- data.frame(type = rep("Reference study", each = 4*B), 
                            dimension = rep(c("beta1", "beta2", "beta3", "beta4"), each = B),
                            states = c(df_pilot[, 1], df_pilot[, 2], df_pilot[, 3], df_pilot[, 4]))
df_exp_long <- data.frame(type = rep("Current study with\n uninformative prior", each = 4*B), 
                          dimension = rep(c("beta1", "beta2", "beta3", "beta4"), each = B),
                          states = c(df_exp[, 1], df_exp[, 2], df_exp[, 3], df_exp[, 4]))
df_combined <- rbind(df_pilot_long, df_exp_long)
df_combined$dimension <- factor(df_combined$dimension,
                                levels =c("beta1", "beta2", "beta3", "beta4"),
                                labels = c(expression(beta[1]), expression(beta[2]),
                                           expression(beta[3]), expression(beta[4])))
custom_colors <- c("Reference study" = "cyan", "Current study with\n uninformative prior" = "blue")
box_plot <- ggplot(df_combined, aes(x = dimension, y = states, fill = type)) +
  geom_boxplot(fatten = 0.8) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]),
                              expression(beta[3]), expression(beta[4]))) +
  xlab("Dimension") +
  ylab("Value") +
  labs(fill = "") +
  scale_fill_manual(values = custom_colors) +
  guides(fill = guide_legend(
    keyheight = unit(4, "line"),  
    label.theme = element_text(size = 25) 
  )) +
  theme(text = element_text(size = 25),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)))
ggsave("/Users/chenyangzhong/Desktop/Data_paper/Figs/app1.pdf", plot = box_plot, device = cairo_pdf, width = 12, height = 7, dpi = 300, bg = "white")

# Functions for graph-enabled MCMC
## Log likelihood function
log_L_func <- function(X, y, n, beta){
  res <- 0
  for(i in 1:n){
    temp <- y[i]*(X[i,]%*%t(beta)) - log(1+exp(X[i,]%*%t(beta)))
    res <- res + temp
  }
  return(res)
}
## Log prior
log_prior <- function(beta){
  res <- log(dmvnorm(beta, mean=rep(0,2), sigma=diag(c(2.5^2,5.3^2),2,2)))
  return(res)
}
## Generalized graph-enabled MCMC
### Parameters: n.chain: number of chains
###             chain.length: length of one chain
###             k, h, hp, rho: algorithm parameters, where
###                hp: standard deviation for Gaussian noise added to \theta_T
###             Theta: prior samples for theta_C
###             theta_init: initial value for theta_T
###             X, y: observation matrix 
mcmc.sim <- function(n.chain=2, chain.length, k, h, hp, rho, Theta, theta_init, X, y){
  
  n <- nrow(X)
  d <- ncol(X)
  B <- nrow(Theta)
  d0 <- ncol(Theta)
  d1 <- d - d0
  
  return_list = NULL
  
  start.init = Sys.time()
  
  dist_vec <- dist(Theta, method = "euclidean", diag = T)
  dist_mat <- as.matrix(dist_vec, B, B)
  diag(dist_mat) <- Inf
  num_idx <- 1:B
  nbh <- num_idx[rank(dist_mat[1,], ties.method="min")<=k]
  graph_edgelist <- cbind(rep(1,length(nbh)), nbh)
  
  for(i in 2:B){
    nbh <- num_idx[rank(dist_mat[i,], ties.method="min")<=k]
    graph_edgelist <- rbind(graph_edgelist, cbind(rep(i,length(nbh)), nbh))
  }
  
  idx <- !duplicated(t(apply(graph_edgelist, 1, sort))) 
  graph_edgelist <- graph_edgelist[idx,]
  
  A <- matrix(0, nrow=B, ncol=B)
  nbh_record <- NULL
  for(i in 1:B){
    nbh_list <- graph_edgelist[graph_edgelist[,1]==i, 2]
    nbh_list <- c(nbh_list, graph_edgelist[graph_edgelist[,2]==i, 1])
    nbh_record <- append(nbh_record, list(nbh_list))
    A[i, nbh_list] <- 1
  }
  DN = rowSums(A)
  
  init.time = Sys.time() - start.init
  
  run.time = rep(0, chain.length)
  
  for(c in 1:n.chain){
    
    start.time <- Sys.time()
    
    alpha_idx <- sample(1:B, 1)
    alpha_tilde <- rmvnorm(1, mean = Theta[alpha_idx,], sigma = diag(h^2,d0))
    alpha_p <- theta_init
    alpha_tilde_vec <- c(alpha_p, alpha_tilde)
    
    run.time[1] = Sys.time() - start.time
    
    t <- 1
    
    while(t < chain.length){
      
      start.time <- Sys.time()
      
      ind <- rbinom(1, 1, rho)
      if(ind == 0){
        beta_idx <- sample(nbh_record[[alpha_idx]], 1)
      } else {
        beta_idx <- sample(1:B, 1)
      }
      eps <- rmvnorm(1, mean = rep(0,d0), sigma=diag(h^2,d0))
      eps_p <- rmvnorm(1, mean=rep(0,d1), sigma=diag(hp^2,d1))
      beta_tilde <- Theta[beta_idx,] + eps
      beta_p <- alpha_p + eps_p
      indicator <- A[alpha_idx, beta_idx]
      temp <- c(log_L_func(X, y, n, t(as.matrix(c(beta_p,beta_tilde)))), 
                log_L_func(X, y, n, t(as.matrix(c(alpha_p,alpha_tilde)))),
                rho/B+indicator*(1-rho)/DN[beta_idx], rho/B+indicator*(1-rho)/DN[alpha_idx], 
                log_prior(beta_p), log_prior(alpha_p))
      
      log_MH_ratio <- log(temp[3])+temp[1]+temp[5]-log(temp[4])-temp[2]-temp[6]
      MH_ratio <- exp(log_MH_ratio)
      
      if(runif(1,0,1) < MH_ratio){
        
        alpha_idx = beta_idx
        alpha_tilde = beta_tilde
        alpha_p = beta_p
        alpha_tilde_vec <- rbind(alpha_tilde_vec, c(beta_p,beta_tilde))
        
      }else{
        alpha_tilde_vec <- rbind(alpha_tilde_vec, c(alpha_p,alpha_tilde))
      }
      t <- t+1
      run.time[t] <- Sys.time()-start.time
    }
    
    res_list <- list("init_time" = init.time, "run_time" = run.time, 
                     "alpha_df" = alpha_tilde_vec)
    return_list[[c]] <- res_list
    print(paste0("Graph-enabled MCMC: Chain ",c))
  }
  return(return_list)
}

# Simulation using graph-enabled MCMC
dic_name <- "/Users/chenyangzhong/Desktop/Data_paper/Application/"
X_dat <- as.matrix(cbind(rep(1,nrow(gallup2)), gallup2[,c("fluid_dummy", "round_num")],
                         gallup2[,"round_num"]*gallup2[,"fluid_dummy"]))
Theta <- Theta_pilot[, 3:4]
num_sample <- dim(Theta)[1]
saveRDS(Theta, paste0(dic_name,"Theta"))
Theta <- readRDS(paste0(dic_name,"Theta"))
n.chain = 3
chain.length = 100000
k <- ceiling(sqrt(B))
h <- 0.01
hp <- 0.05
rho <- 0.5
theta.init <- c(mean(Theta_exp[,1]), mean(Theta_exp[,2]))
result <- mcmc.sim(n.chain = n.chain, chain.length = chain.length, k=k,  h=h, hp=hp, rho=rho, Theta = Theta, theta_init = theta.init, X = X_dat, y = gallup2$decision..0.D.1.C.)
list.save(result, paste0(dic_name, "result.rdata"))

# Visualizing the results from graph-enabled MCMC
dic_name <- "/Users/chenyangzhong/Desktop/Data_paper/Application/"
result <- list.load(paste0(dic_name,"result.rdata"))
samples <- result[[1]]$alpha_df
df_pilot <- data.frame(Theta_pilot[,3:4])
df_exp <- data.frame(Theta_exp[,3:4])
df_sample <- data.frame(samples[50001:100000,3:4])
df_pilot_long <- data.frame(type = rep("Reference study", each = 2*B), 
                            dimension = rep(c("beta3", "beta4"), each = B),
                            states = c(df_pilot[, 1], df_pilot[, 2]))
df_exp_long <- data.frame(type = rep("Current study with\n uninformative prior", each = 2*B), 
                          dimension = rep(c("beta3", "beta4"), each = B),
                          states = c(df_exp[, 1], df_exp[, 2]))
df_sample_long <- data.frame(type = rep("Current study with\n informative prior", each = 2*50000), 
                             dimension = rep(c("beta3", "beta4"), each = 50000),
                             states = c(df_sample[, 1], df_sample[, 2]))
df_combined <- rbind(df_pilot_long, df_exp_long, df_sample_long)
df_combined$dimension <- factor(df_combined$dimension,
                                levels =c("beta3", "beta4"),
                                labels = c(expression(beta[3]), expression(beta[4])))
custom_colors <- c("Reference study" = "cyan", "Current study with\n uninformative prior" = "blue", "Current study with\n informative prior" = "red")
box_plot <- ggplot(df_combined, aes(x = dimension, y = states, fill = type)) +
  geom_boxplot(fatten = 0.8) +
  scale_x_discrete(labels = c(expression(beta[3]), expression(beta[4]))) +
  xlab("Dimension") +
  ylab("Value") +
  labs(fill = "")+
  scale_fill_manual(values = custom_colors) +
  guides(fill = guide_legend(
    keyheight = unit(5, "line"),  
    label.theme = element_text(size = 28) 
  )) +
  theme(text = element_text(size = 28),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)))
ggsave("/Users/chenyangzhong/Desktop/Data_paper/Figs/app2.pdf", plot = box_plot, device = cairo_pdf, width = 10, height = 7.5, dpi = 300, bg = "white")

data1 <- df_pilot
data2 <- df_exp
data3 <- df_sample
data1$group <- "Group 1"
data2$group <- "Group 2"
data3$group <- "Group 3"
colnames(data1)[1:2]=c("x","y")
colnames(data2)[1:2]=c("x","y")
colnames(data3)[1:2]=c("x","y")
color_map <- c("Reference study" = "cyan2", "Current study with\n uninformative prior" = "blue", "Current study with\n informative prior" = "red")
p <- ggplot(data1, aes(x=x,y=y)) +
  coord_equal() +
  geom_point(data = data1, aes(x, y, color = "Reference study"), size = 1) +
  geom_density_2d(data = data2, aes(color = "Current study with\n uninformative prior"), h = 0.1) + 
  geom_density_2d(data = data3, aes(color = "Current study with\n informative prior"), h = 0.1) +
  scale_colour_manual(name = "", values = color_map) +
  guides(color = guide_legend(nrow = 3,
                              keyheight = unit(5, "line"),
                              keywidth = unit(3, "line"),
                              override.aes = list(size = 2),
                              label.theme = element_text(size = 28))) +
  labs(title = "", x = expression(beta[3]), y = expression(beta[4]), color = NULL) +
  theme(text = element_text(size = 28),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        plot.title = element_text(hjust = 0.5))
ggsave("/Users/chenyangzhong/Desktop/Data_paper/Figs/app3.pdf", plot = p, device = cairo_pdf, width = 10, height = 7.5, dpi = 300, bg = "white")



