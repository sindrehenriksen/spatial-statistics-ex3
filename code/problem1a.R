setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(tidyverse)
library(MASS)
library(reshape)
seismic <- as.vector(read.table('../data/seismic.dat',row.names=NULL))
seismic.df <- melt(vapply(seq(1,75),rep, rep(1,75), times= 75))
seismic.df$value <- seismic$V1

ggplot(seismic.df,aes(x= X1,y= X2)) + 
  geom_tile(aes(fill=value))+
  scale_fill_gradient(low ="white",
                      high ="black")+
  theme_minimal()
