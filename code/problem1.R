setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

library(tidyverse)
library(MASS)
library(reshape)
library(ggplot2)

## ---- a
seismic = as.vector(
  read.table('../data/seismic.dat', row.names=NULL))
seismic.df = melt(vapply(seq(1, 75), rep, rep(1,75), times=75))
seismic.df$value = seismic$V1

seismic.obs =
  ggplot(seismic.df, aes(x=X1, y=X2, fill=value)) +
  geom_raster() +
  scale_fill_gradient2() +
  theme_minimal()
print(seismic.obs)
#ggsave("../figures/seismic_data.pdf", plot = seismic.obs,
#       width=5, height=3.5, units="in")

## ---- b
n = nrow(seismic)
p_di_li0 = dnorm(seismic$V1, mean=0.02, sd=0.06)
p_di_li1 = dnorm(seismic$V1, mean=0.08, sd=0.06)
p_di = 0.5 * (p_di_li0 + p_di_li1)
p_li1_di = 0.5 * p_di_li1 / p_di
b_post_reals =
  replicate(10, as.integer(rbernoulli(n, p=p_li1_di)))

# Plot simulations
b_data =
  as_tibble(b_post_reals) %>%
  bind_cols(dplyr::select(seismic.df, X1, X2)) %>%
  gather(sim, value, -X1, -X2)

b_sims_plot =
  ggplot(b_data, aes(x=X1, y=X2, fill=factor(value))) +
  facet_wrap(~sim) +
  geom_raster() +
  scale_fill_brewer() +
  theme_minimal() +
  labs(fill="Value")
ggsave("../figures/b_sims.pdf", plot=b_sims_plot,
       width=6.5, height=5, units="in")

# Plot mean
b_mean_plot =
  dplyr::select(seismic.df, X1, X2) %>%
  add_column(mean=p_li1_di) %>%
  ggplot(aes(x=X1, y=X2, fill=mean)) +
  geom_raster() +
  theme_minimal()
ggsave("../figures/b_mean.pdf", plot=b_mean_plot,
       width=5, height=3.5, units="in")

# Plot variance
b_var_plot =
  dplyr::select(seismic.df, X1, X2) %>%
  add_column(var=p_li1_di * (1 - p_li1_di)) %>%
  ggplot(aes(x=X1, y=X2, fill=var)) +
  geom_raster() +
  theme_minimal()
ggsave("../figures/b_var.pdf", plot=b_var_plot,
       width=5, height=3.5, units="in")

# Plot MMAP
b_mmap_plot =
  dplyr::select(seismic.df, X1, X2) %>%
  add_column(mmap=as.integer(p_li1_di > 0.5)) %>%
  ggplot(aes(x=X1, y=X2, fill=factor(mmap))) +
  geom_raster() +
  scale_fill_brewer() +
  theme_minimal() +
  labs(fill="MMAP")
ggsave("../figures/b_mmap.pdf", plot=b_mmap_plot,
       width=5, height=3.5, units="in")

## ---- c

#Plotting data from comparable domain
complit_mat = as.matrix(read.table('../data/complit.dat', row.names=NULL))
complit_vec = as.vector(t(complit_mat))
complit.df = melt(vapply(seq(1, 66), rep, rep(1,66), times=66))
complit.df$value = complit_vec

complit.obs =
  ggplot(complit.df, aes(x=X1, y=X2, fill=value)) +
  geom_raster() +
  #scale_fill_gradient2() +
  theme_minimal()
print(complit.obs)
ggsave("../figures/complit_data.pdf", plot = complit.obs,
        width=5, height=3.5, units="in")

#Finding neighbors
find_neighbors = function(values, i, dim, wrap = FALSE){
  n = length(values)
  jumps = c(-dim, -1, 1, dim)
  valid_jump = ifelse((i + jumps) < 1 | (i + jumps) > n, FALSE,TRUE)
  
  if((i %% dim) == 1){
    valid_jump[2] <- FALSE
  }
  if((i %% dim) == 0){
    valid_jump[3] <- FALSE
  }
  
  neighbors = i + jumps[valid_jump]
  
  #Add functionality for wrapping boundary conditions
}

#log likelihood for 1 position
loglik1 = function(beta, values, i , dim){
  possible_vals = c(0,1)
  neighb= find_neighbors(values, i, dim, wrap = FALSE)
  norm_const = exp(beta*sum(possible_vals[1]==values[neighb])) + 
    exp(beta*sum(possible_vals[2]==values[neighb]))
  return(beta*sum(values[i]==values[neighb]) - log(norm_const))
}

loglikAll = function(beta, values, dim){
  n = length(values)
  log_p = rep(0,n)
  for(i in 1:n){
    log_p[i] = loglik1(beta, values, i, dim)
  }
  return(sum(log_p))
}

#Estimating beta
beta_start = 1
beta_hat = optim(beta_start, loglikAll, 
                 values = complit_vec, dim = 66, 
                 method = 'Brent', upper = 6, lower = 0,
                 control = list(fnscale = -1))


