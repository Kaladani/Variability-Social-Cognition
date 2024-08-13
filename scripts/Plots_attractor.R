library(ggplot2)
library(viridis)
library(quantreg)
library(ggpubr)

theme_set(theme_bw(base_size = 10,
                   base_family = "sans"))


## Simulate attraction and repulsion 

set.seed(4242)

n <- 200

pre <- rnorm(n, 0, 100)

hist(pre)

treatment <- c(rep(0, n/2), rep(1,n/2))
treatment <- factor(treatment, labels = c(0,1))
treatment

post <- pre + (treatment == 1) * (-pre)  + rnorm(n, 0, 25)



mydata <- data.frame(pre = pre, treatment = treatment, post = post)
p_a <- mydata |> 
  ggplot(aes(x = post, color = treatment, fill = treatment, group = treatment)) + 
  geom_density(alpha=0.4, position="identity",
                 binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T) +
  labs(x="",y="density") +
  scale_y_continuous(breaks = c(0,0.005,0.01,0.015,0.02), limits = c(0,0.02)) +
  scale_x_continuous(breaks = c(-300,-200,-100,0,100,200,300), limits = c(-300,300)) +
  scale_color_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1) +
  scale_fill_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1)


## ideal point somewhere else

n <- 200

pre <- rnorm(n, 0, 100)

hist(pre)

treatment <- c(rep(0, n/2), rep(1,n/2))
treatment <- factor(treatment, labels = c(0,1))
treatment


post <- pre + (treatment == 1) * (200-pre) + rnorm(n, 0, 25)



mydata <- data.frame(pre = pre, treatment = treatment, post = post)
p_b <- mydata |> 
  ggplot(aes(x = post, color = treatment, fill = treatment, group = treatment)) + 
  geom_density(alpha=0.4, position="identity",
               binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T) +
  labs(x="",y="density") +
  scale_y_continuous(breaks = c(0,0.005,0.01,0.015,0.02), limits = c(0,0.02)) +
  scale_x_continuous(breaks = c(-300,-200,-100,0,100,200,300), limits = c(-300,300)) +
  scale_color_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1) +
  scale_fill_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1)


library(quantreg)

quant <- rq(formula = post ~ treatment, tau = c(0.1,0.25,0.5,0.75,0.9), data = mydata)
summary(quant)
anova(quant)

# repulson center



n <- 200

pre <- rnorm(n, 0, 10)
pre <- pre*2

hist(pre)

treatment <- c(rep(0, n/2), rep(1,n/2))
treatment <- factor(treatment, labels = c(0,1))
treatment


post <- pre + (treatment == 1) * 50*sign(pre)/sqrt(abs(pre)) + rnorm(n, 0, 5)
post <- post*2

mydata <- data.frame(pre = pre, treatment = treatment, post = post)
p_c <- mydata |> 
  ggplot(aes(x = post, color = treatment, fill = treatment, group = treatment)) + 
  geom_density(alpha=0.4, position="identity",
               binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T) +
  labs(x="",y="density") +
  scale_y_continuous(breaks = c(0,0.005,0.01,0.015,0.02), limits = c(0,0.02)) +
  scale_x_continuous(breaks = c(-300,-200,-100,0,100,200,300), limits = c(-150,150)*2) +
  scale_color_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1) +
  scale_fill_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1)

p_c

library(quantreg)

quant <- rq(formula = post ~ treatment, tau = c(0.1,0.25,0.5,0.75,0.9), data = mydata)
summary(quant)
anova(quant)

# summary(rq(formula = post ~ pre*treatment, tau = c(0.1,0.25,0.5,0.75,0.9), data = mydata))

## repulsion asymmetric floor


n <- 200

rnormt <- function(n, range, mu, s) {
  
  # range is a vector of two values
  
  F.a <- pnorm(min(range), mean = mu, sd = s)
  F.b <- pnorm(max(range), mean = mu, sd = s)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qnorm(u, mean = mu, sd = s)
  
}

pre <- rnormt(n, c(0, Inf), 50, 50)
pre <- pre

hist(pre)

treatment <- c(rep(0, n/2), rep(1,n/2))
treatment <- factor(treatment, labels = c(0,1))
treatment


post <- pre + (treatment == 1) * (200+rnorm(n/2, 0, 10))/(pre)^(1/3)
post <- post
library(ggplot2)
library(viridis)

theme_set(theme_bw(base_size = 10,
                   base_family = "sans"))

mydata <- data.frame(pre = pre, treatment = treatment, post = post)

tapply(mydata$post, mydata$treatment, summary)


p_d <- mydata |> 
  ggplot(aes(x = post, color = treatment, fill = treatment, group = treatment)) + 
  geom_density(alpha=0.4, position="identity",
               binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T,
               bounds = c(0,Inf)) +
  labs(x="",y="density") +
  scale_y_continuous(breaks = c(0,0.005,0.01,0.015,0.02), limits = c(0,0.02)) +
  scale_x_continuous(breaks = c(-300,-200,-100,0,100,200,300), limits = c(-150,150)*2) +
  scale_color_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1) +
  scale_fill_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1)




quant <- rq(formula = post ~ treatment, tau = c(0.1,0.25,0.5,0.75,0.9), data = mydata)
summary(quant)
anova(quant)


p_attract <- ggarrange(p_a, p_c, p_b, p_d, ncol = 2, nrow =2, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(c)", "(b)", "(d)"), vjust = 1.5,
          font.label = list(size = 8, color = "black", face = "bold", family = "sans"))

p_attract

ggsave(file="plots/Plot_Attractor-Repulsor-Data.tiff", p_attract, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/Plot_Attractor-Repulsor-Data.jpeg", p_attract, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)
