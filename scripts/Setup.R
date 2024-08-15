rm(list=ls())

pck_names <- c("ggplot2", 
               "viridis", 
               "ggpubr", 
               "patchwork", 
               "tidyverse",
               "dglm",
               "quantreg",
               "statmod",
               "rsimsum",
               "MASS",
               "marginaleffects",
               "broom",
               "texreg")


install.packages(pck_names, dep=T)

lapply(pck_names, require, character.only = TRUE)