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


# Parameter values
sigma = 1
tau = 2
n = 10
h = 1
mu1 = c(4,0)
mu2 = c(-4,0)
mu3 = c(0,4)
B = 100
n.chain <- 3
chain.length <- 10000
k <- ceiling(sqrt(B))
rho <- 0.5
sigma_eps <- 0.5

# Read data
dic_name <- "/Users/chenyangzhong/Desktop/Data_paper/GM/"
X <- readRDS(paste0(dic_name,"X.rds"))
Theta0 <- readRDS(paste0(dic_name,"Theta0.rds"))
theta <- readRDS(paste0(dic_name,"theta.rds"))
res_1 <- list.load(paste0(dic_name,"result.rdata"))
res_2 <- list.load(paste0(dic_name,"result_mrw.rdata"))
result_true <- list.load(paste0(dic_name,"result_true.rdata"))
samples_1 <- res_1[[1]]$alpha_df
samples_2 <- res_2[[1]]$alpha_df
samples_true <- result_true[[1]]$alpha_df

# Wasserstein distance computation
set.seed(0)
prior.mean = colMeans(Theta0)
prior.cov = cov(Theta0)
mean_gk = as.vector(solve(solve(prior.cov)+n*tau^(-2)*diag(1,2))%*%(solve(prior.cov)%*%prior.mean+tau^(-2)*colSums(X)))
cov_gk = solve(solve(prior.cov)+n*tau^(-2)*diag(1,2))
W2G = 0
W2M = 0
W2gk = 0 
for(c in 1:n.chain){
  sample_true <- result_true[[c]]$alpha_df
  sample_gk <- mvrnorm(n=chain.length/2, mu=mean_gk, Sigma=cov_gk)
  sample_1 <- res_1[[c]]$alpha_df
  sample_2 <- res_2[[c]]$alpha_df
  W2G = W2G + wasserstein(pp(sample_true[(chain.length/2+1):chain.length,]), pp(sample_1[(chain.length/2+1):chain.length,]), p=2)
  W2M = W2M + wasserstein(pp(sample_true[(chain.length/2+1):chain.length,]), pp(sample_2[(chain.length/2+1):chain.length,]), p=2)
  W2gk = W2gk + wasserstein(pp(sample_true[(chain.length/2+1):chain.length,]), pp(sample_gk), p=2)
}
W2G = W2G/(n.chain)
W2M = W2M/(n.chain)
W2gk = W2gk/(n.chain)
saveRDS(W2G,paste0(dic_name,"W2G.rds"))
saveRDS(W2M,paste0(dic_name,"W2M.rds"))
saveRDS(W2gk,paste0(dic_name,"W2gk.rds"))

# Running time comparison
init.time_1 <- 0
run.time_1 <- 0
run.time_2 <- 0
for(c in 1:n.chain){
  init.time_1 <- init.time_1 + res_1[[c]]$init_time
  run.time_1 <- run.time_1 + res_1[[c]]$run_time
  run.time_2 <- run.time_2 + res_2[[c]]$run_time  
}
init.time_1 = init.time_1/n.chain
run.time_1 = run.time_1/n.chain
run.time_2 = run.time_2/n.chain

temp = init.time_1
tim_1 <- rep(0, chain.length)
for(i in 1:chain.length){
  temp = temp + run.time_1[i]
  tim_1[i] = temp 
}

temp = 0
tim_2 <- rep(0, chain.length)
for(i in 1:chain.length){
  temp = temp + run.time_2[i]
  tim_2[i] = temp 
}

tim_1[chain.length]/chain.length
tim_2[chain.length]/chain.length

# MCMC convergence diagnostics
chain.prep <- function(dat, n.chain){
  return_list = NULL
  for(c in 1:n.chain){
    return_list[[c]] <- as.mcmc(dat[[c]]$alpha_df)
  }
  return(return_list)
}
gmcmc_chain <- chain.prep(res_1, n.chain)
gmcmc_converge_test <- gelman.diag(gmcmc_chain)$mpsrf
print(gmcmc_converge_test)
mrw_chain <- chain.prep(res_2, n.chain)
mrw_converge_test <- gelman.diag(mrw_chain)$mpsrf
print(mrw_converge_test)

## Rhat plot
gp_g = gelman.plot(gmcmc_chain)
gp_m = gelman.plot(mrw_chain)
df_g = data.frame(bind_rows(as.data.frame(gp_g[["shrink"]][,,1]), 
                            as.data.frame(gp_g[["shrink"]][,,2])), 
                  q=rep(dimnames(gp_g[["shrink"]])[[3]], each=nrow(gp_g[["shrink"]][,,1])),
                  last.iter=rep(gp_g[["last.iter"]], length(gp_g)))
colnames(df_g)[1:2] <- c("\u03B8[1]", "\u03B8[2]")
Rhat_g = ggplot(melt(df_g, c("q","last.iter"), value.name="shrink_factor"), 
                aes(last.iter, shrink_factor, colour=q, linetype=q)) + 
  geom_hline(yintercept=1, colour="grey30", lwd=0.2) +
  geom_line(size = 2.5) +
  facet_grid(~variable, labeller = label_parsed) +
  ggtitle("Graph-enabled MCMC") +
  labs(x="Last iteration in chain", y="Shrink factor",
       colour="Quantile", linetype="Quantile") +
  scale_linetype_manual(values=c(2,1)) +
  scale_y_continuous(limits = c(1, 3.2)) +
  theme(text = element_text(size=60), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        legend.spacing.y = unit(1.5, 'cm'),
        legend.text = element_text(size=60)) +
  guides(colour = guide_legend(keywidth = 5, keyheight = 5),  
         linetype = guide_legend(keywidth = 5, keyheight = 5))
df_m = data.frame(bind_rows(as.data.frame(gp_m[["shrink"]][,,1]), 
                            as.data.frame(gp_m[["shrink"]][,,2])), 
                  q=rep(dimnames(gp_m[["shrink"]])[[3]], each=nrow(gp_m[["shrink"]][,,1])),
                  last.iter=rep(gp_m[["last.iter"]], length(gp_m)))
colnames(df_m)[1:2] <- c("\u03B8[1]", "\u03B8[2]")
Rhat_m = ggplot(melt(df_m, c("q","last.iter"), value.name="shrink_factor"), 
                aes(last.iter, shrink_factor, colour=q, linetype=q)) + 
  geom_hline(yintercept=1, colour="grey30", lwd=0.2) +
  geom_line(size = 2.5) +
  facet_grid(~variable, labeller = label_parsed) +
  ggtitle("Metropolis random walk")+
  labs(x="Last iteration in chain", y="Shrink factor",
       colour="Quantile", linetype="Quantile") +
  scale_linetype_manual(values=c(2,1)) +
  scale_y_continuous(limits = c(1, 3.2)) +
  theme(text = element_text(size=60), plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        legend.spacing.y = unit(1.5, 'cm'),
        legend.text = element_text(size=60)) +
  guides(colour = guide_legend(keywidth = 5, keyheight = 5),  
         linetype = guide_legend(keywidth = 5, keyheight = 5))
Rhat <- ggarrange(Rhat_g, Rhat_m, nrow = 2)
ggsave("/Users/chenyangzhong/Desktop/Data_paper/Figs/GMM_Rhat.pdf", plot = Rhat, device = cairo_pdf, width = 32, height = 20, dpi = 300, bg = "white")

## Traceplot
df_g <- data.frame(type = rep("Graph-enabled MCMC", each = 2*chain.length), 
                   steps = rep(1:chain.length, 2),
                   dimension = rep(c("\u03B8[1]", "\u03B8[2]"), each = chain.length),
                   states = c(samples_1[, 1], samples_1[, 2]))
df_m <- data.frame(type = rep("Metropolis random walk", each = 2*chain.length), 
                   steps = rep(1:chain.length, 2),
                   dimension = rep(c("\u03B8[1]", "\u03B8[2]"), each = chain.length),
                   states = c(samples_2[, 1], samples_2[, 2]))
df_combined <- rbind(df_g, df_m)
traceplot = ggplot(df_combined, aes(x = steps, y = states)) + 
  geom_line() +
  facet_grid(vars(dimension), vars(type), scales = "free_y", labeller = labeller(dimension = label_parsed)) +
  xlab("Number of steps") +
  ylab("Value") +
  theme(text = element_text(size=80), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.text = element_text(size=80))
ggsave("/Users/chenyangzhong/Desktop/Data_paper/Figs/traceplot.pdf", plot = traceplot, device = cairo_pdf, width = 35, height = 18, dpi = 600, bg = "white")

## ACF_plot
bacf_g_1 <- acf(samples_1[,1], plot = FALSE)
bacf_df_g_1 <- with(bacf_g_1, data.frame(lag, acf))
bacf_df_g_1 <- bacf_df_g_1 %>%
  mutate(type = rep("Graph-enabled MCMC", nrow(bacf_df_g_1)), dimension = rep("\u03B8[1]", nrow(bacf_df_g_1)))

bacf_g_2 <- acf(samples_1[,2], plot = FALSE)
bacf_df_g_2 <- with(bacf_g_2, data.frame(lag, acf))
bacf_df_g_2 <- bacf_df_g_2 %>%
  mutate(type = rep("Graph-enabled MCMC", nrow(bacf_df_g_2)), dimension = rep("\u03B8[2]", nrow(bacf_df_g_2)))

bacf_m_1 <- acf(samples_2[,1], plot = FALSE)
bacf_df_m_1 <- with(bacf_m_1, data.frame(lag, acf))
bacf_df_m_1 <- bacf_df_m_1 %>%
  mutate(type = rep("Metropolis random walk", nrow(bacf_df_m_1)), dimension = rep("\u03B8[1]", nrow(bacf_df_m_1)))

bacf_m_2 <- acf(samples_2[,2], plot = FALSE)
bacf_df_m_2 <- with(bacf_m_2, data.frame(lag, acf))
bacf_df_m_2 <- bacf_df_m_2 %>%
  mutate(type = rep("Metropolis random walk", nrow(bacf_df_m_2)), dimension = rep("\u03B8[2]", nrow(bacf_df_m_2)))

bacf_df_combined <- rbind(bacf_df_g_1, bacf_df_g_2, bacf_df_m_1, bacf_df_m_2)

ACFplot = ggplot(bacf_df_combined, aes(x = lag, y = acf)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  facet_grid(vars(dimension), vars(type), scales = "free_y", labeller = labeller(dimension = label_parsed)) +
  xlab("Lag") +
  ylab("ACF") +
  theme(text = element_text(size=70), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.text = element_text(size=70))
ggsave("/Users/chenyangzhong/Desktop/Data_paper/Figs/acf_GMM.pdf", plot = ACFplot, device = cairo_pdf, width = 22, height = 12, dpi = 600)

## Effective sample size
ess_g <-rep(0,2)
ess_m <-rep(0,2)
for(c in 1:n.chain){
  ess_g = ess_g + effectiveSize(as.mcmc(res_1[[c]]$alpha_df[(chain.length/2+1):chain.length,]))
  ess_m = ess_m + effectiveSize(as.mcmc(res_2[[c]]$alpha_df[(chain.length/2+1):chain.length,]))
}
ess_g/n.chain
ess_m/n.chain

# Visualizing true and estimated posteriors
## Prior samples and likelihood
x_seq = seq(-8, 8, 0.1)
y_seq = seq(-4, 8, 0.1)

true.prior_func <- function(x, y){
  prob = 1/3*dmvnorm(c(x, y), mean=mu1, sigma=sigma^2*diag(2))+
    1/3*dmvnorm(c(x, y), mean=mu2, sigma=sigma^2*diag(2))+
    1/3*dmvnorm(c(x, y), mean=mu3, sigma=sigma^2*diag(2))
  return(prob)
}
likelihood_func <- function(x, y){
  prob = 1
  for(l in 1:n){
    prob = prob*dmvnorm(X[l,], mean = c(x, y), sigma = tau^2*diag(2))
  }
  return(prob)
}

form_data <- expand_grid(x = x_seq, y = y_seq)

temp <- NULL
for(i in 1:nrow(form_data)){
  temp = c(temp,  likelihood_func(as.numeric(form_data[i,1]), as.numeric(form_data[i,2])))
}
temp = as.matrix(temp, ncol = 1)
### Likelihood grid values
l_data <- mutate(form_data, z=temp)

temp <- NULL
for(i in 1:nrow(form_data)){
  temp = c(temp,  true.prior_func(as.numeric(form_data[i,1]), as.numeric(form_data[i,2])))
}
temp = as.matrix(temp, ncol = 1)
### True prior grid values
true.prior_data <- mutate(form_data, z=temp)
### True parameter
annotation_data <- data.frame(x = theta[1], y = theta[2])

## Posteriors
x_seq_new = seq(2, 6, 0.1)
y_seq_new = seq(-1, 2, 0.1)
prior.mean = colMeans(Theta0)
prior.cov = cov(Theta0)

kde_func <- function(x, y){
  prob = 0
  for(k in 1:B){
    prob = prob + 1/(2*pi)*exp(-((x-Theta0[k,1])^2+(y-Theta0[k,2])^2)/(2*h^2))
  }
  prob = prob/(B*h^2)
  return(prob)
}
gk_func <- function(x,y){
  prob = dmvnorm(c(x, y), mean=prior.mean, sigma=prior.cov)
  return(prob)
}
true.posterior_func <- function(x, y){
  prob = true.prior_func(x,y)*likelihood_func(x,y)
  return(prob)
}

form_data_new <- expand_grid(x = x_seq_new, y = y_seq_new)

temp <- NULL
for(i in 1:nrow(form_data_new)){
  temp = c(temp,  true.posterior_func(as.numeric(form_data_new[i,1]), as.numeric(form_data_new[i,2])))
}
temp = as.matrix(temp, ncol = 1)
### True posterior grid values
true.posterior_data <- mutate(form_data_new, z=temp)

temp <- NULL
for(i in 1:nrow(form_data_new)){
  temp = c(temp,  kde_func(as.numeric(form_data_new[i,1]), as.numeric(form_data_new[i,2]))*likelihood_func(as.numeric(form_data_new[i,1]), as.numeric(form_data_new[i,2])))
}
temp = as.matrix(temp, ncol = 1)
### KDE posterior grid values
kde_data <- mutate(form_data_new, z=temp)

temp <- NULL
for(i in 1:nrow(form_data_new)){
  temp = c(temp,  gk_func(as.numeric(form_data_new[i,1]), as.numeric(form_data_new[i,2]))*likelihood_func(as.numeric(form_data_new[i,1]), as.numeric(form_data_new[i,2])))
}
temp = as.matrix(temp, ncol = 1)
### Gaussian posterior grid values
g_data <- mutate(form_data_new, z=temp)

Theta0_df <- data.frame(Theta0)
colnames(Theta0_df) <- c("x", "y")

color_map <- c("Likelihood" = "blue", "True prior" = "red", "Prior samples" = "black", "True parameter" = "orange")
plot1 <- ggplot(l_data, aes(x, y)) +  
  geom_contour(aes(z = z, color = "Likelihood")) +       
  geom_contour(data = true.prior_data, aes(z = z, color = "True prior"), linetype = "dashed") +
  geom_point(data = Theta0_df, aes(x, y, color = "Prior samples")) +
  geom_point(data = annotation_data, aes(x, y, color = "True parameter"), shape = 17, size = 3.5) +
  scale_color_manual(values = color_map) +
  guides(color = guide_legend(nrow = 2)) +
  labs(title = "Prior samples and likelihood", x = expression(theta[1]), y = expression(theta[2]), color = NULL) +
  theme(text = element_text(size = 26), plot.title = element_text(hjust = 0.5), legend.position = "right",
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        legend.text = element_text(size = 28),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "white", color = NA))

color_map <- c("Kernel density estimate" = "black", "Gaussian estimate"="purple", "True posterior" = "red", "True parameter" = "orange")
plot2 <- ggplot(kde_data, aes(x, y)) + 
  geom_contour(aes(z = z, color = "Gaussian estimate")) +
  geom_contour(aes(z = z, color = "Kernel density estimate")) +       
  geom_contour(data = true.posterior_data, aes(z = z, color = "True posterior"), linetype = "dashed") +
  geom_point(data = annotation_data, aes(x, y, color = "True parameter"), shape = 17, size = 3.5) +
  scale_color_manual(values = color_map) +
  guides(color = guide_legend(nrow = 2)) +
  labs(title = "Kernel density estimate", x = expression(theta[1]), y = expression(theta[2]), color = NULL)+
  theme(text = element_text(size = 26), plot.title = element_text(hjust = 0.5), legend.position = "right",
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        legend.text = element_text(size = 28),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "white", color = NA))
plot3 <- ggplot(g_data, aes(x,y))+
  geom_contour(aes(z=z, color = "Gaussian estimate"))+
  geom_contour(data = true.posterior_data, aes(z=z, color = "True posterior"), linetype = "dashed")+
  geom_point(data = annotation_data, aes(x,y, color = "True parameter"), shape = 17, size = 3.5)+
  scale_color_manual(values = color_map) +
  guides(color = guide_legend(nrow = 2)) +
  labs(title = "Gaussian estimate", x = expression(theta[1]), y = expression(theta[2]), color = NULL)+
  theme(text = element_text(size = 26), plot.title = element_text(hjust = 0.5), legend.position = "right",
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        legend.text = element_text(size = 28),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "white", color = NA))

fig1 <- ggarrange(plot1, ggarrange(plot2, plot3, common.legend = TRUE, legend = "bottom"), nrow = 1, widths = c(1, 2), legend = "bottom")
ggsave("/Users/chenyangzhong/Desktop/Data_paper/Figs/GMM_contour.pdf", plot = fig1, device = cairo_pdf, width = 20, height = 7.5, dpi = 300, bg = "white")

## Comparing posterior distributions
df_g = data.frame(x=samples_1[(chain.length/2+1):chain.length,1], y=samples_1[(chain.length/2+1):chain.length,2])
df_m = data.frame(x=samples_2[(chain.length/2+1):chain.length,1], y=samples_2[(chain.length/2+1):chain.length,2])
df_true = data.frame(x=samples_true[(chain.length/2+1):chain.length,1], y=samples_true[(chain.length/2+1):chain.length,2])
p1 <- ggplot(df_g, aes(x=x,y=y)) + coord_equal()
p1 <- p1 + geom_density_2d(aes(color = "Estimated"), h = 1.15) +
  geom_density_2d(data = df_true, aes(color = "True"), h = 1.15, linetype = "dashed")+
  xlab(expression(theta[1])) +
  ylab(expression(theta[2])) +
  ggtitle("Graph-enabled MCMC") +
  theme(text = element_text(size = 35), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8))) +
  scale_color_manual(name = "Posterior", values = c("Estimated" = "blue", "True" = "red")) +
  guides(color = guide_legend(keyheight = unit(2.5, "line"),
                              keywidth = unit(2.5, "line"),
                              label.theme = element_text(size = 35)))
p2 <- ggplot(df_m, aes(x=x,y=y)) + coord_equal()
p2 <- p2 + geom_density_2d(aes(color = "Estimated"), h = 1)+
  geom_density_2d(data = df_true, aes(color = "True"), h = 1.15, linetype = "dashed")+
  xlab(expression(theta[1])) +
  ylab(expression(theta[2])) +
  ggtitle("Metropolis random walk") +
  theme(text = element_text(size = 35), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8))) +
  scale_color_manual(name = "Posterior", values = c("Estimated" = "blue", "True" = "red")) +
  guides(color = guide_legend(keyheight = unit(2.5, "line"),
                              keywidth = unit(2.5, "line"),
                              label.theme = element_text(size = 35)))
p <- ggarrange(p1, p2, common.legend = TRUE, legend = "right") 

df_g_long <- data.frame(type = rep("Graph-enabled MCMC", each = chain.length), 
                        dimension = rep(c("x", "y"), each = chain.length/2),
                        states = c(samples_1[(chain.length/2+1):chain.length, 1], samples_1[(chain.length/2+1):chain.length, 2]))
df_m_long <- data.frame(type = rep("Metropolis random walk", each = chain.length), 
                        dimension = rep(c("x", "y"), each = chain.length/2),
                        states = c(samples_2[(chain.length/2+1):chain.length, 1], samples_2[(chain.length/2+1):chain.length, 2]))
df_true_long <- data.frame(type = rep("True posterior", each = chain.length), 
                           dimension = rep(c("x", "y"), each = chain.length/2),
                           states = c(samples_true[(chain.length/2+1):chain.length, 1], samples_true[(chain.length/2+1):chain.length, 2]))
df_combined <- rbind(df_g_long, df_m_long, df_true_long)
df_combined$dimension <- factor(df_combined$dimension,
                                levels = c("x", "y"),
                                labels = c(expression(theta[1]), expression(theta[2])))
q <- ggplot(df_combined, aes(x = dimension, y = states, fill = type)) +
  geom_boxplot() +
  scale_x_discrete(labels = c(expression(theta[1]), expression(theta[2]))) +
  xlab("Dimension") +
  ylab("Value") +
  labs(fill = "Algorithm")+
  theme_minimal() +
  theme(text = element_text(size = 35),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8))) +
  guides(fill = guide_legend(keyheight = unit(2.5, "line"),
                              keywidth = unit(2.5, "line"),
                              label.theme = element_text(size = 35)))

Comparison <- ggarrange(p, q, nrow = 2) 
ggsave("/Users/chenyangzhong/Desktop/Data_paper/Figs/Comparing.pdf", plot = Comparison, device = cairo_pdf, width = 20, height = 15, dpi = 300, bg = "white")


