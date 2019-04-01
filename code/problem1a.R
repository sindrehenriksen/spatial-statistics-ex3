setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(tidyverse)
library(MASS)
library(reshape)
seismic <- as.vector(read.table('../data/seismic.dat',row.names=NULL))
seismic.df <- melt(vapply(seq(1,75),rep, rep(1,75), times= 75))
seismic.df$value <- seismic$V1

seismic.obs <- ggplot(seismic.df,aes(x= X1,y= X2)) + 
  geom_tile(aes(fill=value))+
  scale_fill_gradient(low ="white",
                      high ="black")+
  theme_minimal()

ggsave("../figures/seismic_data.pdf", plot = seismic.obs, 
       device = NULL, path = NULL, scale = 1, width = 5.5, 
       height = 5.5, units = "in", dpi = 300, limitsize = TRUE)