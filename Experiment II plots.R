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
library(Cairo)

# Setting parameters
n.chain = 3
chain.length = 10000
B = 10000
n.iter = 40000

# Loading data
## d=2
d = 2
dic_name <- paste0("./OUD/d=", d, ",")
samps_ohio <- readRDS(paste0(dic_name,"samps_ohio.rds"))
Theta_ohio_2 <- samps_ohio[[1]][(n.iter-B+1):n.iter,]
samps_ny <- readRDS(paste0(dic_name,"samps_ny.rds"))
Theta_ny_2 <- samps_ny[[1]][(n.iter-B+1):n.iter,]
samps_true <- readRDS(paste0(dic_name,"samps_true.rds"))
res_2_true <- samps_true
Theta_true_2 <- samps_true[[1]][(n.iter-B+1):n.iter,]
dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
res_2_G <- list.load(paste0(dic_name_B, "result.rdata"))
res_2_M <- list.load(paste0(dic_name_B, "result_mrw.rdata"))
samples_2_G <- res_2_G[[1]]$alpha_df
samples_2_M <- res_2_M[[1]]$alpha_df

## d=6
d = 6
dic_name <- paste0("./OUD/d=", d, ",")
samps_ohio <- readRDS(paste0(dic_name,"samps_ohio.rds"))
Theta_ohio_6 <- samps_ohio[[1]][(n.iter-B+1):n.iter,]
samps_ny <- readRDS(paste0(dic_name,"samps_ny.rds"))
Theta_ny_6 <- samps_ny[[1]][(n.iter-B+1):n.iter,]
samps_true <- readRDS(paste0(dic_name,"samps_true.rds"))
res_6_true <- samps_true
Theta_true_6 <- samps_true[[1]][(n.iter-B+1):n.iter,]
dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
res_6_G <- list.load(paste0(dic_name_B, "result.rdata"))
res_6_M <- list.load(paste0(dic_name_B, "result_mrw.rdata"))
samples_6_G <- res_6_G[[1]]$alpha_df
samples_6_M <- res_6_M[[1]]$alpha_df

## d=10
d = 10
dic_name <- paste0("./OUD/d=", d, ",")
samps_ohio <- readRDS(paste0(dic_name,"samps_ohio.rds"))
Theta_ohio_10 <- samps_ohio[[1]][(n.iter-B+1):n.iter,]
samps_ny <- readRDS(paste0(dic_name,"samps_ny.rds"))
Theta_ny_10 <- samps_ny[[1]][(n.iter-B+1):n.iter,]
samps_true <- readRDS(paste0(dic_name,"samps_true.rds"))
res_10_true <- samps_true
Theta_true_10 <- samps_true[[1]][(n.iter-B+1):n.iter,]
dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
res_10_G <- list.load(paste0(dic_name_B, "result.rdata"))
res_10_M <- list.load(paste0(dic_name_B, "result_mrw.rdata"))
samples_10_G <- res_10_G[[1]]$alpha_df
samples_10_M <- res_10_M[[1]]$alpha_df

# 2-Wasserstein distance computation
## d=2
W2G_2 = 0
W2M_2 = 0
for(c in 1:n.chain){
  sample_true <- res_2_true[[c]]
  sample_1 <- res_2_G[[c]]$alpha_df
  sample_2 <- res_2_M[[c]]$alpha_df
  W2G_2 = W2G_2 + wasserstein(pp(sample_true[(n.iter-chain.length/2+1):n.iter,]), pp(sample_1[(chain.length/2+1):chain.length,]), p=2)
  W2M_2 = W2M_2 + wasserstein(pp(sample_true[(n.iter-chain.length/2+1):n.iter,]), pp(sample_2[(chain.length/2+1):chain.length,]), p=2)
}
W2G_2 = W2G_2/(n.chain)
W2M_2 = W2M_2/(n.chain)
### Graph-enabled MCMC
print(W2G_2)
### Metropolis random walk
print(W2M_2)
## d=6
W2G_6 = 0
W2M_6 = 0
for(c in 1:n.chain){
  sample_true <- res_6_true[[c]]
  sample_1 <- res_6_G[[c]]$alpha_df
  sample_2 <- res_6_M[[c]]$alpha_df
  W2G_6 = W2G_6 + wasserstein(pp(sample_true[(n.iter-chain.length/2+1):n.iter,]), pp(sample_1[(chain.length/2+1):chain.length,]), p=2)
  W2M_6 = W2M_6 + wasserstein(pp(sample_true[(n.iter-chain.length/2+1):n.iter,]), pp(sample_2[(chain.length/2+1):chain.length,]), p=2)
}
W2G_6 = W2G_6/(n.chain)
W2M_6 = W2M_6/(n.chain)
### Graph-enabled MCMC
print(W2G_6)
### Metropolis random walk
print(W2M_6)
## d=10
W2G_10 = 0
W2M_10 = 0
for(c in 1:n.chain){
  sample_true <- res_10_true[[c]]
  sample_1 <- res_10_G[[c]]$alpha_df
  sample_2 <- res_10_M[[c]]$alpha_df
  W2G_10 = W2G_10 + wasserstein(pp(sample_true[(n.iter-chain.length/2+1):n.iter,]), pp(sample_1[(chain.length/2+1):chain.length,]), p=2)
  W2M_10 = W2M_10 + wasserstein(pp(sample_true[(n.iter-chain.length/2+1):n.iter,]), pp(sample_2[(chain.length/2+1):chain.length,]), p=2)
}
W2G_10 = W2G_10/(n.chain)
W2M_10 = W2M_10/(n.chain)
### Graph-enabled MCMC
print(W2G_10)
### Metropolis random walk
print(W2M_10)

# Comparing distributions
color_map <- c("Ohio posterior samples" = "cyan2", "New York posterior with uninformative prior" = "red", 
               "New York posterior with informative prior via graph-enabled MCMC" = "black")
## d=2
df_ohio_2 = data.frame(x=Theta_ohio_2[,1], y=Theta_ohio_2[,2])
df_ny_2 = data.frame(x=Theta_ny_2[,1], y=Theta_ny_2[,2])
df_G_2 = data.frame(x=samples_2_G[(chain.length/2+1):chain.length,1], y=samples_2_G[(chain.length/2+1):chain.length,2])
P_2 <- ggplot(df_ohio_2, aes(x=x,y=y)) + coord_equal() +
  geom_point(aes(x, y, color = "Ohio posterior samples"), size = 1) +
  geom_density_2d(data = df_ny_2, aes(color = "New York posterior with uninformative prior"), h = 0.2, linetype = "dashed") + 
  geom_density_2d(data = df_G_2, aes(color = "New York posterior with informative prior via graph-enabled MCMC"), h = 0.2) +
  scale_color_manual(name = "", values = color_map) +
  guides(color = guide_legend(nrow = 3,
                              keyheight = unit(3, "line"),
                              keywidth = unit(3, "line"))) +
  labs(title = "d=2", x = expression(beta[1]), y = expression(beta[2]), color = NULL) +
  theme(text = element_text(size = 30), 
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),
        legend.position = "bottom", legend.text = element_text(size = 28))
## d=6
df_ohio_6 = data.frame(x=Theta_ohio_6[,1], y=Theta_ohio_6[,2])
df_ny_6 = data.frame(x=Theta_ny_6[,1], y=Theta_ny_6[,2])
df_G_6 = data.frame(x=samples_6_G[(chain.length/2+1):chain.length,1], y=samples_6_G[(chain.length/2+1):chain.length,2])
P_6 <- ggplot(df_ohio_6, aes(x=x,y=y)) + coord_equal() +
  geom_point(aes(x, y, color = "Ohio posterior samples"), size = 1) +
  geom_density_2d(data = df_ny_6, aes(color = "New York posterior with uninformative prior"), h = 0.2, linetype = "dashed") + 
  geom_density_2d(data = df_G_6, aes(color = "New York posterior with informative prior via graph-enabled MCMC"), h = 0.2) +
  scale_color_manual(name = "", values = color_map) +
  guides(color = guide_legend(nrow = 3,
                              keyheight = unit(3, "line"),
                              keywidth = unit(3, "line"))) +
  labs(title = "d=6", x = expression(beta[1]), y = expression(beta[2]), color = NULL) +
  theme(text = element_text(size = 30),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 10)), 
        legend.position = "bottom", legend.text = element_text(size = 28))
## d=10
df_ohio_10 = data.frame(x=Theta_ohio_10[,1], y=Theta_ohio_10[,2])
df_ny_10 = data.frame(x=Theta_ny_10[,1], y=Theta_ny_10[,2])
df_G_10 = data.frame(x=samples_10_G[(chain.length/2+1):chain.length,1], y=samples_10_G[(chain.length/2+1):chain.length,2])
P_10 <- ggplot(df_ohio_10, aes(x=x,y=y)) + coord_equal() +
  geom_point(aes(x, y, color = "Ohio posterior samples"), size = 1) +
  geom_density_2d(data = df_ny_10, aes(color = "New York posterior with uninformative prior"), h = 0.25, linetype = "dashed") + 
  geom_density_2d(data = df_G_10, aes(color = "New York posterior with informative prior via graph-enabled MCMC"), h = 0.25) +
  scale_color_manual(name = "", values = color_map) +
  guides(color = guide_legend(nrow = 3,
                              keyheight = unit(3, "line"),
                              keywidth = unit(3, "line"))) +
  labs(title = "d=10", x = expression(beta[1]), y = expression(beta[2]), color = NULL) +
  theme(text = element_text(size = 30), 
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),
        legend.position = "bottom", legend.text = element_text(size = 28))
## Combining the figures
P <- ggarrange(P_2, P_6, P_10, common.legend = TRUE, legend = "bottom", ncol=3, nrow=1, align = "hv")
ggsave("./Figs/OUD_contour.pdf", plot = P, device = cairo_pdf, width = 20, height = 7.5, dpi = 300, bg = "white")

# Running time comparison
## d=2
d = 2
dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
res_2_10000_G <- list.load(paste0(dic_name_B, "result.rdata"))
res_2_10000_M <- list.load(paste0(dic_name_B, "result_mrw.rdata"))
init.time_2_10000_G <- 0
run.time_2_10000_G <- 0
run.time_2_10000_M <- 0
for(c in 1:n.chain){
  init.time_2_10000_G <- init.time_2_10000_G + res_2_10000_G[[c]]$init_time
  run.time_2_10000_G <- run.time_2_10000_G + res_2_10000_G[[c]]$run_time
  run.time_2_10000_M <- run.time_2_10000_M + res_2_10000_M[[c]]$run_time  
}
init.time_2_10000_G = init.time_2_10000_G/n.chain
run.time_2_10000_G = run.time_2_10000_G/n.chain
run.time_2_10000_M = run.time_2_10000_M/n.chain
ave.time_2_10000_G = as.numeric((init.time_2_10000_G+sum(run.time_2_10000_G))/chain.length)
ave.time_2_10000_M = as.numeric(sum(run.time_2_10000_M)/chain.length)
print(ave.time_2_10000_G)
print(ave.time_2_10000_M)
## d=6
d = 6
dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
res_6_10000_G <- list.load(paste0(dic_name_B, "result.rdata"))
res_6_10000_M <- list.load(paste0(dic_name_B, "result_mrw.rdata"))
init.time_6_10000_G <- 0
run.time_6_10000_G <- 0
run.time_6_10000_M <- 0
for(c in 1:n.chain){
  init.time_6_10000_G <- init.time_6_10000_G + res_6_10000_G[[c]]$init_time
  run.time_6_10000_G <- run.time_6_10000_G + res_6_10000_G[[c]]$run_time
  run.time_6_10000_M <- run.time_6_10000_M + res_6_10000_M[[c]]$run_time  
}
init.time_6_10000_G = init.time_6_10000_G/n.chain
run.time_6_10000_G = run.time_6_10000_G/n.chain
run.time_6_10000_M = run.time_6_10000_M/n.chain
ave.time_6_10000_G = as.numeric((init.time_6_10000_G+sum(run.time_6_10000_G))/chain.length)
ave.time_6_10000_M = as.numeric(sum(run.time_6_10000_M)/chain.length)
print(ave.time_6_10000_G)
print(ave.time_6_10000_M)
## d=10
d = 10
dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
res_10_10000_G <- list.load(paste0(dic_name_B, "result.rdata"))
res_10_10000_M <- list.load(paste0(dic_name_B, "result_mrw.rdata"))
init.time_10_10000_G <- 0
run.time_10_10000_G <- 0
run.time_10_10000_M <- 0
for(c in 1:n.chain){
  init.time_10_10000_G <- init.time_10_10000_G + res_10_10000_G[[c]]$init_time
  run.time_10_10000_G <- run.time_10_10000_G + res_10_10000_G[[c]]$run_time
  run.time_10_10000_M <- run.time_10_10000_M + res_10_10000_M[[c]]$run_time  
}
init.time_10_10000_G = init.time_10_10000_G/n.chain
run.time_10_10000_G = run.time_10_10000_G/n.chain
run.time_10_10000_M = run.time_10_10000_M/n.chain
ave.time_10_10000_G = as.numeric((init.time_10_10000_G+sum(run.time_10_10000_G))/chain.length)
ave.time_10_10000_M = as.numeric(sum(run.time_10_10000_M)/chain.length)
print(ave.time_10_10000_G)
print(ave.time_10_10000_M)

# ACF plots
## d=2
### Graph-enabled MCMC, beta1
bacf_2_G_1 <- acf(samples_2_G[,1], plot = FALSE)
bacf_df_2_G_1 <- with(bacf_2_G_1, data.frame(lag, acf))
bacf_df_2_G_1 <- bacf_df_2_G_1 %>%
  mutate(type = rep("Graph-enabled MCMC", nrow(bacf_df_2_G_1)), dimension = rep("\u03B2[1]", nrow(bacf_df_2_G_1)))
### Graph-enabled MCMC, beta2
bacf_2_G_2 <- acf(samples_2_G[,2], plot = FALSE)
bacf_df_2_G_2 <- with(bacf_2_G_2, data.frame(lag, acf))
bacf_df_2_G_2 <- bacf_df_2_G_2 %>%
  mutate(type = rep("Graph-enabled MCMC", nrow(bacf_df_2_G_2)), dimension = rep("\u03B2[2]", nrow(bacf_df_2_G_2)))
### Metropolis random walk, beta1
bacf_2_M_1 <- acf(samples_2_M[,1], plot = FALSE)
bacf_df_2_M_1 <- with(bacf_2_M_1, data.frame(lag, acf))
bacf_df_2_M_1 <- bacf_df_2_M_1 %>%
  mutate(type = rep("Metropolis random walk", nrow(bacf_df_2_M_1)), dimension = rep("\u03B2[1]", nrow(bacf_df_2_M_1)))
### Metropolis random walk, beta2
bacf_2_M_2 <- acf(samples_2_M[,2], plot = FALSE)
bacf_df_2_M_2 <- with(bacf_2_M_2, data.frame(lag, acf))
bacf_df_2_M_2 <- bacf_df_2_M_2 %>%
  mutate(type = rep("Metropolis random walk", nrow(bacf_df_2_M_2)), dimension = rep("\u03B2[2]", nrow(bacf_df_2_M_2)))
### plot
bacf_df_2_combined <- rbind(bacf_df_2_G_1, bacf_df_2_G_2, bacf_df_2_M_1, bacf_df_2_M_2)
ACFplot_2 = ggplot(bacf_df_2_combined, aes(x = lag, y = acf)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  facet_grid(vars(dimension), vars(type), scales = "free_y", labeller = labeller(dimension = label_parsed)) +
  xlab("Lag") +
  ylab("ACF") +
  theme(text = element_text(size=70), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.text = element_text(size=70))
ggsave("./Figs/acf d=2.pdf", plot = ACFplot_2, device = cairo_pdf, width = 22, height = 12, dpi = 600)

## d=6
### Graph-enabled MCMC, beta1
bacf_6_G_1 <- acf(samples_6_G[,1], plot = FALSE)
bacf_df_6_G_1 <- with(bacf_6_G_1, data.frame(lag, acf))
bacf_df_6_G_1 <- bacf_df_6_G_1 %>%
  mutate(type = rep("Graph-enabled MCMC", nrow(bacf_df_6_G_1)), dimension = rep("\u03B2[1]", nrow(bacf_df_6_G_1)))
### Graph-enabled MCMC, beta2
bacf_6_G_2 <- acf(samples_6_G[,2], plot = FALSE)
bacf_df_6_G_2 <- with(bacf_6_G_2, data.frame(lag, acf))
bacf_df_6_G_2 <- bacf_df_6_G_2 %>%
  mutate(type = rep("Graph-enabled MCMC", nrow(bacf_df_6_G_2)), dimension = rep("\u03B2[2]", nrow(bacf_df_6_G_2)))
### Metropolis random walk, beta1
bacf_6_M_1 <- acf(samples_6_M[,1], plot = FALSE)
bacf_df_6_M_1 <- with(bacf_6_M_1, data.frame(lag, acf))
bacf_df_6_M_1 <- bacf_df_6_M_1 %>%
  mutate(type = rep("Metropolis random walk", nrow(bacf_df_6_M_1)), dimension = rep("\u03B2[1]", nrow(bacf_df_6_M_1)))
### Metropolis random walk, beta2
bacf_6_M_2 <- acf(samples_6_M[,2], plot = FALSE)
bacf_df_6_M_2 <- with(bacf_6_M_2, data.frame(lag, acf))
bacf_df_6_M_2 <- bacf_df_6_M_2 %>%
  mutate(type = rep("Metropolis random walk", nrow(bacf_df_6_M_2)), dimension = rep("\u03B2[2]", nrow(bacf_df_6_M_2)))
### plot
bacf_df_6_combined <- rbind(bacf_df_6_G_1, bacf_df_6_G_2, bacf_df_6_M_1, bacf_df_6_M_2)
ACFplot_6 = ggplot(bacf_df_6_combined, aes(x = lag, y = acf)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  facet_grid(vars(dimension), vars(type), scales = "free_y", labeller = labeller(dimension = label_parsed)) +
  xlab("Lag") +
  ylab("ACF") +
  theme(text = element_text(size=70), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.text = element_text(size=70))
ggsave("./Figs/acf d=6.pdf", plot = ACFplot_6, device = cairo_pdf, width = 22, height = 12, dpi = 600)

## d=10
### Graph-enabled MCMC, beta1
bacf_10_G_1 <- acf(samples_10_G[,1], plot = FALSE)
bacf_df_10_G_1 <- with(bacf_10_G_1, data.frame(lag, acf))
bacf_df_10_G_1 <- bacf_df_10_G_1 %>%
  mutate(type = rep("Graph-enabled MCMC", nrow(bacf_df_10_G_1)), dimension = rep("\u03B2[1]", nrow(bacf_df_10_G_1)))
### Graph-enabled MCMC, beta2
bacf_10_G_2 <- acf(samples_10_G[,2], plot = FALSE)
bacf_df_10_G_2 <- with(bacf_10_G_2, data.frame(lag, acf))
bacf_df_10_G_2 <- bacf_df_10_G_2 %>%
  mutate(type = rep("Graph-enabled MCMC", nrow(bacf_df_10_G_2)), dimension = rep("\u03B2[2]", nrow(bacf_df_10_G_2)))
### Metropolis random walk, beta1
bacf_10_M_1 <- acf(samples_10_M[,1], plot = FALSE)
bacf_df_10_M_1 <- with(bacf_10_M_1, data.frame(lag, acf))
bacf_df_10_M_1 <- bacf_df_10_M_1 %>%
  mutate(type = rep("Metropolis random walk", nrow(bacf_df_10_M_1)), dimension = rep("\u03B2[1]", nrow(bacf_df_10_M_1)))
### Metropolis random walk, beta2
bacf_10_M_2 <- acf(samples_10_M[,2], plot = FALSE)
bacf_df_10_M_2 <- with(bacf_10_M_2, data.frame(lag, acf))
bacf_df_10_M_2 <- bacf_df_10_M_2 %>%
  mutate(type = rep("Metropolis random walk", nrow(bacf_df_10_M_2)), dimension = rep("\u03B2[2]", nrow(bacf_df_10_M_2)))
### plot
bacf_df_10_combined <- rbind(bacf_df_10_G_1, bacf_df_10_G_2, bacf_df_10_M_1, bacf_df_10_M_2)
ACFplot_10 = ggplot(bacf_df_10_combined, aes(x = lag, y = acf)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  facet_grid(vars(dimension), vars(type), scales = "free_y", labeller = labeller(dimension = label_parsed)) +
  xlab("Lag") +
  ylab("ACF") +
  theme(text = element_text(size=70), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.text = element_text(size=70))
ggsave("./Figs/acf d=10.pdf", plot = ACFplot_10, device = cairo_pdf, width = 22, height = 12, dpi = 600)

# MCMC convergence diagnostics
chain.prep <- function(dat, n.chain){
  return_list = NULL
  for(c in 1:n.chain){
    return_list[[c]] <- as.mcmc(dat[[c]]$alpha_df)
  }
  return(return_list)
}
## MPSRF, d=2
chain_2_G <- chain.prep(res_2_G, n.chain)
converge_test_2_G <- gelman.diag(chain_2_G)$mpsrf
print(converge_test_2_G)
chain_2_M <- chain.prep(res_2_M, n.chain)
converge_test_2_M <- gelman.diag(chain_2_M)$mpsrf
print(converge_test_2_M)
## MPSRF, d=6
chain_6_G <- chain.prep(res_6_G, n.chain)
converge_test_6_G <- gelman.diag(chain_6_G)$mpsrf
print(converge_test_6_G)
chain_6_M <- chain.prep(res_6_M, n.chain)
converge_test_6_M <- gelman.diag(chain_6_M)$mpsrf
print(converge_test_6_M)
## MPSRF, d=10
chain_10_G <- chain.prep(res_10_G, n.chain)
converge_test_10_G <- gelman.diag(chain_10_G)$mpsrf
print(converge_test_10_G)
chain_10_M <- chain.prep(res_10_M, n.chain)
converge_test_10_M <- gelman.diag(chain_10_M)$mpsrf
print(converge_test_10_M)

# Rhat plot
chain.prep2 <- function(dat, n.chain){
  return_list = NULL
  for(c in 1:n.chain){
    return_list[[c]] <- as.mcmc(dat[[c]]$alpha_df[,1:2])
  }
  return(return_list)
}
chain_2_G <- chain.prep2(res_2_G, n.chain)
chain_2_M <- chain.prep2(res_2_M, n.chain)
chain_6_G <- chain.prep2(res_6_G, n.chain)
chain_6_M <- chain.prep2(res_6_M, n.chain)
chain_10_G <- chain.prep2(res_10_G, n.chain)
chain_10_M <- chain.prep2(res_10_M, n.chain)
## d=2
gp_2_G = gelman.plot(chain_2_G)
gp_2_M = gelman.plot(chain_2_M)
df_2_G = data.frame(bind_rows(as.data.frame(gp_2_G[["shrink"]][,,1]), 
                            as.data.frame(gp_2_G[["shrink"]][,,2])), 
                  q=rep(dimnames(gp_2_G[["shrink"]])[[3]], each=nrow(gp_2_G[["shrink"]][,,1])),
                  last.iter=rep(gp_2_G[["last.iter"]], length(gp_2_G)))
colnames(df_2_G)[1:2] <- c("\u03B2[1]", "\u03B2[2]")
Rhat_2_G = ggplot(melt(df_2_G, c("q","last.iter"), value.name="shrink_factor"), 
                aes(last.iter, shrink_factor, colour=q, linetype=q)) + 
  geom_hline(yintercept=1, colour="grey30", lwd=0.2) +
  geom_line(size = 1.5) +
  facet_grid(~variable, labeller = label_parsed) +
  ggtitle("Graph-enabled MCMC") +
  labs(x="Last Iteration in Chain", y="Shrink Factor",
       colour="Quantile", linetype="Quantile") +
  scale_linetype_manual(values=c(2,1)) +
  theme(text = element_text(size=40), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.spacing.y = unit(1.5, 'cm'),
        legend.text = element_text(size=40)) +
  guides(colour = guide_legend(keywidth = 3, keyheight = 3),  
         linetype = guide_legend(keywidth = 3, keyheight = 3))
df_2_M = data.frame(bind_rows(as.data.frame(gp_2_M[["shrink"]][,,1]), 
                            as.data.frame(gp_2_M[["shrink"]][,,2])), 
                  q=rep(dimnames(gp_2_M[["shrink"]])[[3]], each=nrow(gp_2_M[["shrink"]][,,1])),
                  last.iter=rep(gp_2_M[["last.iter"]], length(gp_2_M)))
colnames(df_2_M)[1:2] <- c("\u03B2[1]", "\u03B2[2]")
Rhat_2_M = ggplot(melt(df_2_M, c("q","last.iter"), value.name="shrink_factor"), 
                aes(last.iter, shrink_factor, colour=q, linetype=q)) + 
  geom_hline(yintercept=1, colour="grey30", lwd=0.2) +
  geom_line(size = 1.5) +
  facet_grid(~variable, labeller = label_parsed) +
  ggtitle("Metropolis random walk")+
  labs(x="Last Iteration in Chain", y="Shrink Factor",
       colour="Quantile", linetype="Quantile") +
  scale_linetype_manual(values=c(2,1)) +
  theme(text = element_text(size=40), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.spacing.y = unit(1.5, 'cm'),
        legend.text = element_text(size=40)) +
  guides(colour = guide_legend(keywidth = 3, keyheight = 3),  
         linetype = guide_legend(keywidth = 3, keyheight = 3))
Rhat_2 <- ggarrange(Rhat_2_G, Rhat_2_M, nrow = 2)
ggsave("./Figs/R d=2.pdf", plot = Rhat_2, device = cairo_pdf, width = 22, height = 12, dpi = 300)
## d=6
gp_6_G = gelman.plot(chain_6_G)
gp_6_M = gelman.plot(chain_6_M)
df_6_G = data.frame(bind_rows(as.data.frame(gp_6_G[["shrink"]][,,1]), 
                              as.data.frame(gp_6_G[["shrink"]][,,2])), 
                    q=rep(dimnames(gp_6_G[["shrink"]])[[3]], each=nrow(gp_6_G[["shrink"]][,,1])),
                    last.iter=rep(gp_6_G[["last.iter"]], length(gp_6_G)))
colnames(df_6_G)[1:2] <- c("\u03B2[1]", "\u03B2[2]")
Rhat_6_G = ggplot(melt(df_6_G, c("q","last.iter"), value.name="shrink_factor"), 
                  aes(last.iter, shrink_factor, colour=q, linetype=q)) + 
  geom_hline(yintercept=1, colour="grey30", lwd=0.2) +
  geom_line(size = 1.5) +
  facet_grid(~variable, labeller = label_parsed) +
  ggtitle("Graph-enabled MCMC") +
  labs(x="Last Iteration in Chain", y="Shrink Factor",
       colour="Quantile", linetype="Quantile") +
  scale_linetype_manual(values=c(2,1)) +
  theme(text = element_text(size=40), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.spacing.y = unit(1.5, 'cm'),
        legend.text = element_text(size=40)) +
  guides(colour = guide_legend(keywidth = 3, keyheight = 3),  
         linetype = guide_legend(keywidth = 3, keyheight = 3))
df_6_M = data.frame(bind_rows(as.data.frame(gp_6_M[["shrink"]][,,1]), 
                              as.data.frame(gp_6_M[["shrink"]][,,2])), 
                    q=rep(dimnames(gp_6_M[["shrink"]])[[3]], each=nrow(gp_6_M[["shrink"]][,,1])),
                    last.iter=rep(gp_6_M[["last.iter"]], length(gp_6_M)))
colnames(df_6_M)[1:2] <- c("\u03B2[1]", "\u03B2[2]")
Rhat_6_M = ggplot(melt(df_6_M, c("q","last.iter"), value.name="shrink_factor"), 
                  aes(last.iter, shrink_factor, colour=q, linetype=q)) + 
  geom_hline(yintercept=1, colour="grey30", lwd=0.2) +
  geom_line(size = 1.5) +
  facet_grid(~variable, labeller = label_parsed) +
  ggtitle("Metropolis random walk")+
  labs(x="Last Iteration in Chain", y="Shrink Factor",
       colour="Quantile", linetype="Quantile") +
  scale_linetype_manual(values=c(2,1)) +
  theme(text = element_text(size=40), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.spacing.y = unit(1.5, 'cm'),
        legend.text = element_text(size=40)) +
  guides(colour = guide_legend(keywidth = 3, keyheight = 3),  
         linetype = guide_legend(keywidth = 3, keyheight = 3))
Rhat_6 <- ggarrange(Rhat_6_G, Rhat_6_M, nrow = 2)
ggsave("./Figs/R d=6.pdf", plot = Rhat_6, device = cairo_pdf, width = 22, height = 12, dpi = 300)
## d=10
gp_10_G = gelman.plot(chain_10_G)
gp_10_M = gelman.plot(chain_10_M)
df_10_G = data.frame(bind_rows(as.data.frame(gp_10_G[["shrink"]][,,1]), 
                              as.data.frame(gp_10_G[["shrink"]][,,2])), 
                    q=rep(dimnames(gp_10_G[["shrink"]])[[3]], each=nrow(gp_10_G[["shrink"]][,,1])),
                    last.iter=rep(gp_10_G[["last.iter"]], length(gp_10_G)))
colnames(df_10_G)[1:2] <- c("\u03B2[1]", "\u03B2[2]")
Rhat_10_G = ggplot(melt(df_10_G, c("q","last.iter"), value.name="shrink_factor"), 
                  aes(last.iter, shrink_factor, colour=q, linetype=q)) + 
  geom_hline(yintercept=1, colour="grey30", lwd=0.2) +
  geom_line(size = 1.5) +
  facet_grid(~variable, labeller = label_parsed) +
  ggtitle("Graph-enabled MCMC") +
  labs(x="Last Iteration in Chain", y="Shrink Factor",
       colour="Quantile", linetype="Quantile") +
  scale_linetype_manual(values=c(2,1))+
  theme(text = element_text(size=40), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.spacing.y = unit(1.5, 'cm'),
        legend.text = element_text(size=40)) +
  guides(colour = guide_legend(keywidth = 3, keyheight = 3),  
         linetype = guide_legend(keywidth = 3, keyheight = 3))
df_10_M = data.frame(bind_rows(as.data.frame(gp_10_M[["shrink"]][,,1]), 
                              as.data.frame(gp_10_M[["shrink"]][,,2])), 
                    q=rep(dimnames(gp_10_M[["shrink"]])[[3]], each=nrow(gp_10_M[["shrink"]][,,1])),
                    last.iter=rep(gp_10_M[["last.iter"]], length(gp_10_M)))
colnames(df_10_M)[1:2] <- c("\u03B2[1]", "\u03B2[2]")
Rhat_10_M = ggplot(melt(df_10_M, c("q","last.iter"), value.name="shrink_factor"), 
                  aes(last.iter, shrink_factor, colour=q, linetype=q)) + 
  geom_hline(yintercept=1, colour="grey30", lwd=0.2) +
  geom_line(size = 1.5) +
  facet_grid(~variable, labeller = label_parsed) +
  ggtitle("Metropolis random walk")+
  labs(x="Last Iteration in Chain", y="Shrink Factor",
       colour="Quantile", linetype="Quantile") +
  scale_linetype_manual(values=c(2,1)) +
  theme(text = element_text(size=40), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.spacing.y = unit(1.5, 'cm'),
        legend.text = element_text(size=40)) +
  guides(colour = guide_legend(keywidth = 3, keyheight = 3),  
         linetype = guide_legend(keywidth = 3, keyheight = 3))
Rhat_10 <- ggarrange(Rhat_10_G, Rhat_10_M, nrow = 2)
ggsave("./Figs/R d=10.pdf", plot = Rhat_10, device = cairo_pdf, width = 22, height = 12, dpi = 300)

# Traceplot
## d=2
df_2_G <- data.frame(type = rep("Graph-enabled MCMC", each = 2*chain.length), 
                   steps = rep(1:chain.length, 2),
                   dimension = rep(c("\u03B2[1]", "\u03B2[2]"), each = chain.length),
                   states = c(samples_2_G[, 1], samples_2_G[, 2]))
df_2_M <- data.frame(type = rep("Metropolis random walk", each = 2*chain.length), 
                   steps = rep(1:chain.length, 2),
                   dimension = rep(c("\u03B2[1]", "\u03B2[2]"), each = chain.length),
                   states = c(samples_2_M[, 1], samples_2_M[, 2]))
df_combined_2 <- rbind(df_2_G, df_2_M)
traceplot_2 = ggplot(df_combined_2, aes(x = steps, y = states)) + 
  geom_line() +
  facet_grid(vars(dimension), vars(type), scales = "free_y", labeller = labeller(dimension = label_parsed)) +
  xlab("Number of steps") +
  ylab("Value") +
  theme(text = element_text(size=70), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.text = element_text(size=70))
ggsave("./Figs/traceplot d=2.pdf", plot = traceplot_2, device = cairo_pdf, width = 30, height = 15, dpi = 600, bg = "white")
## d=6
df_6_G <- data.frame(type = rep("Graph-enabled MCMC", each = 2*chain.length), 
                     steps = rep(1:chain.length, 2),
                     dimension = rep(c("\u03B2[1]", "\u03B2[2]"), each = chain.length),
                     states = c(samples_6_G[, 1], samples_6_G[, 2]))
df_6_M <- data.frame(type = rep("Metropolis random walk", each = 2*chain.length), 
                     steps = rep(1:chain.length, 2),
                     dimension = rep(c("\u03B2[1]", "\u03B2[2]"), each = chain.length),
                     states = c(samples_6_M[, 1], samples_6_M[, 2]))
df_combined_6 <- rbind(df_6_G, df_6_M)
traceplot_6 = ggplot(df_combined_6, aes(x = steps, y = states)) + 
  geom_line() +
  facet_grid(vars(dimension), vars(type), scales = "free_y", labeller = labeller(dimension = label_parsed)) +
  xlab("Number of steps") +
  ylab("Value") +
  theme(text = element_text(size=70), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.text = element_text(size=70))
ggsave("./Figs/traceplot d=6.pdf", plot = traceplot_6, device = cairo_pdf, width = 30, height = 15, dpi = 600, bg = "white")
## d=10
df_10_G <- data.frame(type = rep("Graph-enabled MCMC", each = 2*chain.length), 
                     steps = rep(1:chain.length, 2),
                     dimension = rep(c("\u03B2[1]", "\u03B2[2]"), each = chain.length),
                     states = c(samples_10_G[, 1], samples_10_G[, 2]))
df_10_M <- data.frame(type = rep("Metropolis random walk", each = 2*chain.length), 
                     steps = rep(1:chain.length, 2),
                     dimension = rep(c("\u03B2[1]", "\u03B2[2]"), each = chain.length),
                     states = c(samples_10_M[, 1], samples_10_M[, 2]))
df_combined_10 <- rbind(df_10_G, df_10_M)
traceplot_10 = ggplot(df_combined_10, aes(x = steps, y = states)) + 
  geom_line() +
  facet_grid(vars(dimension), vars(type), scales = "free_y", labeller = labeller(dimension = label_parsed)) +
  xlab("Number of steps") +
  ylab("Value") +
  theme(text = element_text(size=70), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.text = element_text(size=70))
ggsave("./Figs/traceplot d=10.pdf", plot = traceplot_10, device = cairo_pdf, width = 30, height = 15, dpi = 600, bg = "white")

# Effective sample size 
## d = 2
ess_2_G <-rep(0,2)
ess_2_M <-rep(0,2)
for(c in 1:n.chain){
  ess_2_G = ess_2_G + effectiveSize(as.mcmc(res_2_G[[c]]$alpha_df[(chain.length/2+1):chain.length,1:2]))
  ess_2_M = ess_2_M + effectiveSize(as.mcmc(res_2_M[[c]]$alpha_df[(chain.length/2+1):chain.length,1:2]))
}
ess_2_G/n.chain
ess_2_M/n.chain
## d = 6
ess_6_G <-rep(0,2)
ess_6_M <-rep(0,2)
for(c in 1:n.chain){
  ess_6_G = ess_6_G + effectiveSize(as.mcmc(res_6_G[[c]]$alpha_df[(chain.length/2+1):chain.length,1:2]))
  ess_6_M = ess_6_M + effectiveSize(as.mcmc(res_6_M[[c]]$alpha_df[(chain.length/2+1):chain.length,1:2]))
}
ess_6_G/n.chain
ess_6_M/n.chain
## d = 10
ess_10_G <-rep(0,2)
ess_10_M <-rep(0,2)
for(c in 1:n.chain){
  ess_10_G = ess_10_G + effectiveSize(as.mcmc(res_10_G[[c]]$alpha_df[(chain.length/2+1):chain.length,1:2]))
  ess_10_M = ess_10_M + effectiveSize(as.mcmc(res_10_M[[c]]$alpha_df[(chain.length/2+1):chain.length,1:2]))
}
ess_10_G/n.chain
ess_10_M/n.chain

# Running time for d=6 and other B values
d = 6
##B=1000
B <- 1000
dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
res_6_1000_G <- list.load(paste0(dic_name_B, "result.rdata"))
res_6_1000_M <- list.load(paste0(dic_name_B, "result_mrw.rdata"))
init.time_6_1000_G <- 0
run.time_6_1000_G <- 0
run.time_6_1000_M <- 0
for(c in 1:n.chain){
  init.time_6_1000_G <- init.time_6_1000_G + res_6_1000_G[[c]]$init_time
  run.time_6_1000_G <- run.time_6_1000_G + res_6_1000_G[[c]]$run_time
  run.time_6_1000_M <- run.time_6_1000_M + res_6_1000_M[[c]]$run_time  
}
init.time_6_1000_G = init.time_6_1000_G/n.chain
run.time_6_1000_G = run.time_6_1000_G/n.chain
run.time_6_1000_M = run.time_6_1000_M/n.chain
ave.time_6_1000_G = as.numeric((init.time_6_1000_G+sum(run.time_6_1000_G))/chain.length)
ave.time_6_1000_M = as.numeric(sum(run.time_6_1000_M)/chain.length)
##B=2500
B <- 2500
dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
res_6_2500_G <- list.load(paste0(dic_name_B, "result.rdata"))
res_6_2500_M <- list.load(paste0(dic_name_B, "result_mrw.rdata"))
init.time_6_2500_G <- 0
run.time_6_2500_G <- 0
run.time_6_2500_M <- 0
for(c in 1:n.chain){
  init.time_6_2500_G <- init.time_6_2500_G + res_6_2500_G[[c]]$init_time
  run.time_6_2500_G <- run.time_6_2500_G + res_6_2500_G[[c]]$run_time
  run.time_6_2500_M <- run.time_6_2500_M + res_6_2500_M[[c]]$run_time  
}
init.time_6_2500_G = init.time_6_2500_G/n.chain
run.time_6_2500_G = run.time_6_2500_G/n.chain
run.time_6_2500_M = run.time_6_2500_M/n.chain
ave.time_6_2500_G = as.numeric((init.time_6_2500_G+sum(run.time_6_2500_G))/chain.length)
ave.time_6_2500_M = as.numeric(sum(run.time_6_2500_M)/chain.length)
##B=5000
B <- 5000
dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
res_6_5000_G <- list.load(paste0(dic_name_B, "result.rdata"))
res_6_5000_M <- list.load(paste0(dic_name_B, "result_mrw.rdata"))
init.time_6_5000_G <- 0
run.time_6_5000_G <- 0
run.time_6_5000_M <- 0
for(c in 1:n.chain){
  init.time_6_5000_G <- init.time_6_5000_G + res_6_5000_G[[c]]$init_time
  run.time_6_5000_G <- run.time_6_5000_G + res_6_5000_G[[c]]$run_time
  run.time_6_5000_M <- run.time_6_5000_M + res_6_5000_M[[c]]$run_time  
}
init.time_6_5000_G = init.time_6_5000_G/n.chain
run.time_6_5000_G = run.time_6_5000_G/n.chain
run.time_6_5000_M = run.time_6_5000_M/n.chain
ave.time_6_5000_G = as.numeric((init.time_6_5000_G+sum(run.time_6_5000_G))/chain.length)
ave.time_6_5000_M = as.numeric(sum(run.time_6_5000_M)/chain.length)
##B=15000
B <- 15000
dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
res_6_15000_G <- list.load(paste0(dic_name_B, "result.rdata"))
res_6_15000_M <- list.load(paste0(dic_name_B, "result_mrw.rdata"))
init.time_6_15000_G <- 0
run.time_6_15000_G <- 0
run.time_6_15000_M <- 0
for(c in 1:n.chain){
  init.time_6_15000_G <- init.time_6_15000_G + res_6_15000_G[[c]]$init_time
  run.time_6_15000_G <- run.time_6_15000_G + res_6_15000_G[[c]]$run_time
  run.time_6_15000_M <- run.time_6_15000_M + res_6_15000_M[[c]]$run_time  
}
init.time_6_15000_G = init.time_6_15000_G/n.chain
run.time_6_15000_G = run.time_6_15000_G/n.chain
run.time_6_15000_M = run.time_6_15000_M/n.chain
ave.time_6_15000_G = as.numeric((init.time_6_15000_G+sum(run.time_6_15000_G))/chain.length)
ave.time_6_15000_M = as.numeric(sum(run.time_6_15000_M)/chain.length)
##B=20000
B <- 20000
dic_name_B <- paste0("./OUD/d=", d, ",B=", B, ",")
res_6_20000_G <- list.load(paste0(dic_name_B, "result.rdata"))
res_6_20000_M <- list.load(paste0(dic_name_B, "result_mrw.rdata"))
init.time_6_20000_G <- 0
run.time_6_20000_G <- 0
run.time_6_20000_M <- 0
for(c in 1:n.chain){
  init.time_6_20000_G <- init.time_6_20000_G + res_6_20000_G[[c]]$init_time
  run.time_6_20000_G <- run.time_6_20000_G + res_6_20000_G[[c]]$run_time
  run.time_6_20000_M <- run.time_6_20000_M + res_6_20000_M[[c]]$run_time  
}
init.time_6_20000_G = init.time_6_20000_G/n.chain
run.time_6_20000_G = run.time_6_20000_G/n.chain
run.time_6_20000_M = run.time_6_20000_M/n.chain
ave.time_6_20000_G = as.numeric((init.time_6_20000_G+sum(run.time_6_20000_G))/chain.length)
ave.time_6_20000_M = as.numeric(sum(run.time_6_20000_M)/chain.length)
## Plotting the running time
dat <- data.frame(
  x = rep(c(1000, 2500, 5000, 10000, 15000, 20000), 2),
  y = c(ave.time_6_1000_G, ave.time_6_2500_G, ave.time_6_5000_G, ave.time_6_10000_G, ave.time_6_15000_G, ave.time_6_20000_G,
        ave.time_6_1000_M, ave.time_6_2500_M, ave.time_6_5000_M, ave.time_6_10000_M, ave.time_6_15000_M, ave.time_6_20000_M),
  Algorithm = rep(c("Graph-enabled MCMC", "Metropolis random walk"), each = 6)
)
timplots <- ggplot(dat, aes(x = x, y = y, color = Algorithm)) +
  geom_line(linewidth = 2) +
  geom_point(size = 6) +
  labs(title = "", x = "B", y="Time per iteration (s)") +
  theme(text = element_text(size=60), legend.spacing.y = unit(1, 'cm')) +
  guides(color = guide_legend(keywidth = 5, keyheight = 5))
ggsave("./Figs/Time.pdf", plot = timplots, device = cairo_pdf, height = 12, width = 20, dpi = 300, bg = "white")


