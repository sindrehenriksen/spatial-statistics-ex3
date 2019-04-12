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
  ggplot(complit.df, aes(x=X1, y=X2, fill=as.factor(value))) +
  geom_raster() +
  scale_fill_brewer() +
  theme_minimal()+
  labs(fill="Data")
print(complit.obs)
ggsave("../figures/complit_data.pdf", plot = complit.obs,
        width=5, height=3.5, units="in")

#Finding neighbors
find_neighbors = function(values, i, dim){
  n = length(values)
  jumps = c(-dim, -1, 1, dim)
  valid_jump = ifelse((i + jumps) < 1 | (i + jumps) > n, FALSE,TRUE)
  
  if((i %% dim) == 1){
    valid_jump[2] = FALSE
  }
  if((i %% dim) == 0){
    valid_jump[3] = FALSE
  }
  
  neighbors = i + jumps[valid_jump]
}

#Neighbors using wrapping boundary conditions
find_neighbors_wrap = function(values, i, dim){
  n = length(values)
  #c(under, left, right, over)
  jumps = c(-dim, -1, 1, dim)
  
  #under
  if((i - dim) < 1){
    jumps[1]  = n-dim
  }
  #over
  if((i + dim) > n){
    jumps[4] = -(n-dim)
  }
  
  #left
  if((i %% dim) == 1){
    jumps[2] = dim-1 
  }
  #right
  if((i %% dim) == 0){
    jumps[3] = -(dim-1)
  }
  
  neighbors = i + jumps
}


#log likelihood for 1 position
loglik1 = function(beta, values, i , dim){
  possible_vals = c(0,1)
  neighb= find_neighbors(values, i, dim)
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

##############################################
#MCMC
##############################################
#Neighbors using wrapping boundary conditions
find_neighbors_wrap = function(values, i, dim){
  n = length(values)
  #c(under, left, right, over)
  jumps = c(-dim, -1, 1, dim)
  
  #under
  if((i - dim) < 1){
    jumps[1]  = n-dim
  }
  #over
  if((i + dim) > n){
    jumps[4] = -(n-dim)
  }
  
  #left
  if((i %% dim) == 1){
    jumps[2] = dim-1 
  }
  #right
  if((i %% dim) == 0){
    jumps[3] = -(dim-1)
  }
  
  neighbors = i + jumps
}

GibbsOneIt = function(dval, lval, dim, beta, sigma_d, mean_d){
  num_change = 0
  n = length(lval)
  #Draw n positions to change (we will have to do a lot of iterations if not)
  loc = ceiling(n*runif(n))

  for(i in loc){
    #Calculate full conditional probability
    
    #Calculate p(di|li)
    p_di = dnorm(dval[i], mean = mean_d, sd = sigma_d)
    
    #calculatep(li|neighbors)
    neighbors = find_neighbors_wrap(lval, i, dim)
    p_li0_upper = exp(beta * sum(lval[neighbors] == 0))
    p_li1_upper = exp(beta * sum(lval[neighbors] == 1))
  
    #Do I have to normalize now?
    p_li0 = p_li0_upper#/(p_li0_upper + p_li1_upper)
    p_li1 = p_li1_upper#/(p_li0_upper + p_li1_upper)
  
    #Calculate p(li|di, neighbors)
    p_lid0_upper = p_li0*p_di[1]
    p_lid1_upper = p_li1*p_di[2]
  
    final_p0 = p_lid0_upper/(p_lid0_upper + p_lid1_upper)
    
    prev_i = lval[i]
    
    #Propose new value
    u = runif(1)
    if(u<final_p0){
      lval[i] = 0
    }else{
      lval[i] = 1
    }
    if(lval[i] != prev_i){
      num_change = num_change + 1
    }
    #The value is always accepted
  }
  
  return(list(l = lval, change = num_change))
  
}

#Initiate values
set.seed(101)

#from a)
seismic = as.vector(
  read.table('../data/seismic.dat', row.names=NULL))
seismic.df = melt(vapply(seq(1, 75), rep, rep(1,75), times=75))
seismic.df$value = seismic$V1

dval = seismic$V1
dim = 75
sigma_d = 0.06
mean_d = c(0.02,0.08)
beta = 1.72
n = length(dval)
l_start = round(runif(n))

#Run MCMC
num_it = 2000
l_MCMC = matrix(0, ncol = num_it, nrow = n)
num_change = rep(0,num_it)
for (b in 1:num_it){
  res = GibbsOneIt(dval, l_start, dim, beta, sigma_d, mean_d)
  num_change[b] = res$change
  l_MCMC[,b] = res$l
  l_start = res$l
}

#Fraction of 0 and 1
#Used to demonstrate convergence
frac_1 = colSums(l_MCMC)/n
frac_0 = 1-frac_1

ggplot(data.frame(frac = c(frac_1,frac_0), id = as.factor(rep(c(1,0), each = num_it)), 
                  x = rep(1:num_it, 2)))+
  geom_path(aes(x = x, y = frac, col = id))+
  labs(col = 'value', x = 'Iteration', y = 'Fraction')
ggsave("../figures/fractionMCMC.pdf", plot=last_plot(),
       width=5, height=3.5, units="in")
  
#Define burn-in period
burnin = 1:500

#Realizations
r_index = as.integer(seq(600, 2000, length.out = 10))
realizations_MCMC = l_MCMC[,r_index]
MCMC_data =
  as_tibble(realizations_MCMC) %>%
  bind_cols(dplyr::select(seismic.df, X1, X2)) %>%
  gather(sim, value, -X1, -X2)


realiz_MCMC_plot =
  ggplot(MCMC_data, aes(x=X1, y=X2, fill=factor(value))) +
  facet_wrap(~sim) +
  geom_raster() +
  scale_fill_brewer() +
  theme_minimal() +
  labs(fill="Value")
ggsave("../figures/realizations_MCMC.pdf", plot=realiz_MCMC_plot,
       width=6.5, height=5, units="in")


X1 = MCMC_data$X1
X2 = MCMC_data$X2
#Posterior mean
PostMean = rowMeans(l_MCMC[,-burnin])
post.mean =
  ggplot(data.frame(Mean = PostMean, X1, X2), aes(x=X1, y=X2, fill=Mean)) +
  geom_raster() +
  theme_minimal()
print(post.mean)
ggsave("../figures/c_mean.pdf", plot=post.mean,
       width=5, height=3.5, units="in")


#Posterior variance
PostVar = apply(l_MCMC[,-burnin], 1, var)
post.var =
  ggplot(data.frame(Variance = PostVar, X1, X2), aes(x=X1, y=X2, fill=Variance)) +
  geom_raster() +
  theme_minimal()
print(post.var)
ggsave("../figures/c_variance.pdf", plot=post.var,
       width=5, height=3.5, units="in")

#MMAP
position_frac1 = rowSums(l_MCMC[,-burnin])/(num_it-length(burnin))
MMAP = as.numeric(position_frac1>0.5)
post.MMAP =
  ggplot(data.frame(MMAP = MMAP, X1, X2), aes(x=X1, y=X2, fill=as.factor(MMAP))) +
  geom_raster() +
  scale_fill_brewer() +
  theme_minimal() +
  labs(fill="MMAP")
print(post.MMAP)
ggsave("../figures/c_MMAP.pdf", plot=post.MMAP,
       width=5, height=3.5, units="in")

