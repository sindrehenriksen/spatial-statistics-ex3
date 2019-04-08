setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

library(tidyverse)
library(MASS)
library(reshape)

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
ggsave("../figures/seismic_data.pdf", plot = seismic.obs,
       width=5.5, height=4, units="in")

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

b_mean_plot =
  dplyr::select(seismic.df, X1, X2) %>%
  add_column(mean=p_li1_di) %>%
  ggplot(aes(x=X1, y=X2, fill=mean)) +
  geom_raster() +
  theme_minimal()
ggsave("../figures/b_mean.pdf", plot=b_mean_plot,
       width=5.5, height=4, units="in")

b_var_plot =
  dplyr::select(seismic.df, X1, X2) %>%
  add_column(var=p_li1_di * (1 - p_li1_di)) %>%
  ggplot(aes(x=X1, y=X2, fill=var)) +
  geom_raster() +
  theme_minimal()
ggsave("../figures/b_var.pdf", plot=b_var_plot,
       width=5.5, height=4, units="in")
