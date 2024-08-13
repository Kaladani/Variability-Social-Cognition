library(tidyverse)
library(ggplot2)
library(viridis)
library(ggpubr)

theme_set(theme_bw(base_size = 10,
                   base_family = "sans"))



noise <- rnorm(1000, 0, 0.1)

assoc_strength <- rnorm(1000, 2, 1.3)

hist(assoc_strength)
summary(assoc_strength)


# -measurement noise array
# -zugrundeliegende "echte" verteilungen (1 x normal, 1 x binomial mit vielen ziehungen)
# -transformation der zugrundeliegenden verteilungen in beobachtbare effekte 
# (normal+noise, logarithmische transformation von binomialergebnis
  
observed_effect_a <- assoc_strength
hist(observed_effect_a)
summary(observed_effect_a)
var(observed_effect_a)

episodes <- rbinom(1000, 100, 0.2)
hist(episodes)
summary(episodes)


observed_effect_mem <- noise*1.5*episodes^(1/2) + 1/2*episodes^(1/2)
hist(observed_effect_mem)
summary(observed_effect_mem)
var(observed_effect_mem)


mydata <- data.frame(association = assoc_strength, memory = episodes) |> 
  pivot_longer(everything(), names_to = "mechanism", values_to = "value")

hist.background <- ggplot(data=mydata, aes(x=value, color=mechanism, fill=mechanism)) +
  geom_histogram(aes(x=value), alpha=0.4, color="black", show.legend = FALSE, 
                 binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T) +
  labs(y="freq", x = "") +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno")

hist.background.final <- hist.background +
  facet_grid(.~mechanism, labeller = as_labeller(
    c(association = "association strength", 
      memory = "memory episodes")),
             scales = "free_x", switch = 'x') +
  theme(strip.background = element_blank(), 
        strip.placement = "outside",
        axis.title.x = element_blank())

hist.background.final

mydata.obs <- data.frame(association = observed_effect_a, memory = observed_effect_mem) |> 
  pivot_longer(everything(), names_to = "mechanism", values_to = "value")

hist.observed <- ggplot(data=mydata.obs, aes(x=value, color=mechanism, fill=mechanism)) +
  geom_histogram(aes(x=value), alpha=0.4, color="black", show.legend = FALSE, 
                 binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T) + #Freedman-Diaconis Rule for binwidth) + 
  labs(y="freq", x = "observed effect") +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno")

hist.observed.final <- hist.observed +
  facet_grid(.~mechanism, labeller = as_labeller(
    c(association = "obs. association effect", 
      memory = "obs. memory effect")),
     switch = 'x') + 
  theme(strip.background = element_blank(), 
        strip.placement = "outside",
        axis.title.x = element_blank())

hist.observed.final

p_truncated <- ggarrange(hist.background.final, hist.observed.final,nrow = 2)

p_truncated

ggsave(file="plots/Plot_Truncated-Data.tiff", p_truncated, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/Plot_Truncated-Data.jpeg", p_truncated, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)
