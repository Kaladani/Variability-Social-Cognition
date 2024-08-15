library(ggplot2)
library(viridis)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(dglm)
library(quantreg)

theme_set(theme_bw(base_size = 10,
                   base_family = "sans"))


rm(list =ls())

rnormt <- function(n, range, mu, s) {
  
  # range is a vector of two values
  
  F.a <- pnorm(min(range), mean = mu, sd = s)
  F.b <- pnorm(max(range), mean = mu, sd = s)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qnorm(u, mean = mu, sd = s)
  
}

rnormc <- function(n, range, mu, s) {
  
  # range is a vector of two values
  
  dist <- rnorm(n, mu, s)
  dist[dist > 1] <- 1
  dist[dist < 0] <- 0
  
  dist  
}

set.seed(13)

# 1 Prototype theory ----
## 1a necessary and sufficient ----

n_obs <- 1000


X <- rnormt(n_obs, c(0,1), 0.5, 0.2)

hist(X)

Y_exp <- 0.5*X + 0.25

Y <- rnormc(n_obs, c(0,1), Y_exp, 0.1)

### Plots ----


plot(X,Y)

mydata_a <- data.frame(X = X, X = Y)
p_a <- ggplot(data = mydata_a, aes(x = X, y = Y)) + 
  geom_point(size=1, alpha=0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  # geom_density_2d(alpha = 0.5) +
  # geom_boxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.3) +
  ylim(c(0,1)) + xlim(c(0,1)) +
  labs(x = "Prototype Valence", y = "Prosocialness") +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_a

ggsave(file="Plot_Prototype-Theory-a-Simulation.tiff", p_a, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Prototype-Theory-a-Simulation.jpeg", p_a, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

Resid2 <- (Y - (Y_exp))^2

plot(X, Resid2)

mydata_a <- data.frame(X = X, Resid2 = Resid2)
p_a_resid2 <- ggplot(data = mydata_a, aes(x = X, y = Resid2)) + 
  geom_point(size=1, alpha=0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  # geom_density_2d(alpha = 0.5) +
  # geom_boxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.3) +
  ylim(c(0,0.15)) + xlim(c(0,1)) +
  labs(x = "Prototype Valence", y = expression(epsilon^2)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_a_resid2

ggsave(file="Plot_Prototype-Theory-a-Simulation-Residuals.tiff", p_a_resid2, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Prototype-Theory-a-Simulation-Residuals.jpeg", p_a_resid2, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

## 1b Attraction ----

n_obs <- 1000

X <- rnormt(n_obs, c(0,1), 0.5, 0.2)

pre <- rnormc(n_obs, c(0,1), 0.25, 0.2)

Y_attractor <- 0.5*X+0.25

Y <- pre + X * (Y_attractor - pre)



Y_exp <- 0.25 + 0.5*X^2


plot(X,Y)


### Plots ----


mydata_b <- data.frame(X = X, Y = Y)

p_b <- ggplot(data = mydata_b, aes(x = X, y = Y)) + 
  geom_point(size=1, alpha=0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  # geom_density_2d(alpha = 0.5) +
  # geom_boxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.3) +
  labs(x = "Prototype Valence", y = "Prosocialness") +
  ylim(c(0,1)) + xlim(c(0,1)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_b

ggsave(file="Plot_Prototype-Theory-b-Simulation.tiff", p_b, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Prototype-Theory-b-Simulation.jpeg", p_b, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


Resid2 <- (Y - Y_exp)^2

plot(X, Resid2)

mydata_b <- data.frame(pre = pre, X = X, Resid2 = Resid2)
p_b_resid2 <- ggplot(data = mydata_b, aes(x = X, y = Resid2)) + 
  geom_point(size=1, alpha=0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  # geom_density_2d(alpha = 0.5) +
  # geom_boxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.3) +
  ylim(c(0,0.15)) + xlim(c(0,1)) +
  labs(x = "Prototype Valence", y = expression(epsilon^2)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_b_resid2

ggsave(file="Plot_Prototype-Theory-b-Simulation-Residuals.tiff", p_b_resid2, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Prototype-Theory-b-Simulation-Residuals.jpeg", p_b_resid2, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

### Test analysis with predictions ----

mydata <- data.frame(X = X, Y = Y)
quantiles <- c(0.1,0.25,0.5,0.75,0.9)
quant <- rq(formula = Y ~ X+I(X^2), tau = quantiles, data = mydata)
quant


#### quantile coef plot ----
# OLS
lm <- lm(data=mydata, 
         formula =  Y ~  X+I(X^2))

ols <- as.data.frame(coef(lm))
ols.ci <- as.data.frame(confint(lm))
ols2 <- cbind(ols, ols.ci)
ols2 <- tibble::rownames_to_column(ols2, var="term") |> 
  dplyr::filter(!grepl("(Intercept)", term))

coef_labels <- list(
  'I(X^2)'=expression(X^2),
  'X'="X"
)

coef_labeller <- function(variable,value){
  return(coef_labels[value])
}

quantile_coef_plot <- quant |> broom::tidy(se.type = "rank", conf.int = TRUE, conf.level = 0.95) |> 
  dplyr::filter(!grepl("(Intercept)", term)) %>%
  ggplot(aes(x=tau,y=estimate))+
    # quantilie results
  geom_point(color="#27408b", size = 1)+ 
  geom_line(color="#27408b")+ 
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.25, fill="#27408b")+
  
  # OLS results
  geom_hline(data = ols2, aes(yintercept= `coef(lm)`), lty=1, color="red")+
  geom_hline(data = ols2, aes(yintercept= `2.5 %`), lty=2, color="red")+
  geom_hline(data = ols2, aes(yintercept= `97.5 %`), lty=2, color="red")+
  scale_x_continuous(breaks=quantiles) +
  ylab("coef") +
  facet_wrap(~term, scales="free", ncol=2, labeller=coef_labeller)

quantile_coef_plot

ggsave(file="Plot_Prototype-Theory-b-Quantile-coef-plot.tiff", quantile_coef_plot, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Prototype-Theory-b-Quantile-coef-plot.jpeg", quantile_coef_plot, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

summary(quant)
anova(quant, joint = FALSE)


vfr <- dglm(data=mydata,
            formula = Y~X+I(X^2),
            dformula = Y~X+I(X^2), 
            family = gaussian(link = "identity"), dlink = "log", method="reml")
summary(vfr)
anova(vfr)


### non-linear quantiles

mydata <- data.frame(X = X, Y = Y)
p_b_quantile <- ggplot(data = mydata, aes(x = X, y = Y)) + 
  geom_point(alpha=0.3, size=1) +
  # geom_boxplot(aes(group= cut_width(X,0.1)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.3) +
  geom_quantile(
    method = "rqss", quantiles = c(0.1,0.25,0.5,0.75,0.9)
    # formula = y ~ x+I(abs(x)*x), quantiles = c(0.1,0.25,0.5,0.75,0.9), 
    , color="darkblue") +
  # geom_density_2d(alpha = 0.5) +
  # geom_boxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.3) +
  labs(x = "Prototype Valence", y = "Prosocialness") +
  ylim(c(0,1)) + xlim(c(0,1)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_b_quantile

ggsave(file="Plot_Prototype-Theory-b-Simulation-Quantile-Plot.tiff", p_b_quantile, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Prototype-Theory-b-Simulation-Quantile-Plot.jpeg", p_b_quantile, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

out50 <- rq(Y ~ X + I(X^2), 
            data=mydata,tau=c(0.5))
out90 <- rq(Y ~ X + I(X^2), 
            data=mydata,tau=c(0.9))
out75 <- rq(Y ~ X + I(X^2), 
            data=mydata,tau=c(0.75))
out25 <- rq(Y ~ X + I(X^2), 
            data=mydata,tau=c(0.25))
out10 <- rq(Y ~ X + I(X^2), 
            data=mydata,tau=c(0.1))

predictdata <- mydata

predictdata[,c("q50_yhat","q50_lo","q50_hi")] <- predict(out50, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("q90_yhat","q90_lo","q90_hi")] <- predict(out90, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("q75_yhat","q75_lo","q75_hi")] <- predict(out75, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("q25_yhat","q25_lo","q25_hi")] <- predict(out25, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("q10_yhat","q10_lo","q10_hi")] <- predict(out10, newdata=predictdata, interval="confidence", level = 0.95)

quant.predict.final <- predictdata |> 
  tidyr::pivot_longer(starts_with("q"),
               names_to = c("tau", ".value"),
               names_pattern = "^(.*)_(.*)$") |> 
  ggplot(aes(x=X, y=Y, color = tau, fill = tau, group = tau)) + 
  # geom_line(aes(x=X, y=yhatq90)) +
  # geom_ribbon(aes(ymin=yhatq90.lo, ymax=yhatq90.hi), alpha=0.3, linetype="dotted") +
  # geom_line(aes(x=X, y=yhatq75)) +
  # geom_ribbon(aes(ymin=yhatq75.lo, ymax=yhatq75.hi), alpha=0.3, linetype="dotted") +
  # geom_line(aes(x=X, y=yhatq50)) +
  # geom_ribbon(aes(ymin=yhatq50.lo, ymax=yhatq50.hi), alpha=0.3, linetype="dotted") +
  # geom_line(aes(x=X, y=yhatq25)) +
  # geom_ribbon(aes(ymin=yhatq25.lo, ymax=yhatq25.hi), alpha=0.3, linetype="dotted") +
  geom_line(aes(x=X, y=yhat)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.3, linetype="dotted") +
  ylim(c(0,1)) + xlim(c(0,1)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno",
                      name = "Quantiles",
                      labels = quantiles) +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno",
                     name = "Quantiles",
                     labels = quantiles) +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

quant.predict.final

ggsave(file="Plot_Prototype-Theory-b-Simulation-Predicted-Quantile-Plot.tiff", quant.predict.final, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Prototype-Theory-b-Simulation-Predicted-Quantile-Plot.jpeg", quant.predict.final, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

predictdglm <- mydata
predictdglm$yhatv <- predict(vfr, dispersion =1, newdata=predictdglm)
predictdglm$upperhatv <- predictdglm$yhatv+qnorm(0.75,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
predictdglm$highesthatv <- predictdglm$yhatv+qnorm(0.9,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
predictdglm$lowerhatv <- predictdglm$yhatv+qnorm(0.25,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
predictdglm$lowesthatv <- predictdglm$yhatv+qnorm(0.1,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))

logvarpredict <- predict(vfr$dispersion.fit, dispersion=2, newdata=predictdglm, se.fit=TRUE)

predictdglm[,"sdhat"] <- sqrt(exp(logvarpredict$fit))
predictdglm[,"sdlo"] <- sqrt(exp(logvarpredict$fit - qnorm(p=0.975,0,1)*logvarpredict$se.fit))
predictdglm[,"sdhi"] <- sqrt(exp(logvarpredict$fit + qnorm(p=0.975,0,1)*logvarpredict$se.fit))

sd.predict.final <- ggplot(predictdglm, aes(x=X, y=sdhat)) +
  geom_line(aes(x=X, y=sdhat)) +
  geom_ribbon(aes(ymin=sdlo, ymax=sdhi), alpha=0.3, linetype="dotted") +
  labs(y=bquote(widehat(sigma)(Y)))+
  xlim(c(0,1)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

sd.predict.final

ggsave(file="Plot_Prototype-Theory-b-Simulation-Predicted-SD-Plot.tiff", sd.predict.final, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Prototype-Theory-b-Simulation-Predicted-SD-Plot.jpeg", sd.predict.final, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

# 2 Exemplar theory ----

## Variable episodes ----

n_obs <- 1000

X_mean <- rnormt(n_obs, c(0,1), 0.5, 0.2)

Z <- runif(n_obs, 0, 1)

X_pers <- sapply(1:n_obs, function(i) {
  rnormc(100, c(0,1), X_mean[i], sqrt(exp(log(0.005)+8*Z[i])))
}
)

X_sample_mean <- apply(X_pers, 2, function(events) {
  sampled_events <- sample(events, 10)
  mean(sampled_events)
}
)

Y_exp <- 0.25 + 0.5 * X_sample_mean

Y <- rnormc(n_obs, c(0,1), Y_exp, 0.05)

# sd(X_sample_mean)
# sd(X_mean)
# hist(X_sample_mean)
# hist(X_mean)
# 
# # X_sample_mean <- rtruncnorm(1000, a=0, b = 1, X_mean, sqrt(Z/10))
# 
# plot(X_mean, X_sample_mean)
# 
# 
# 
# plot(Z, X_sample_mean)


# 
# plot(X_sample_mean, Y)
# plot(Z, Y)
# 
# plot(X_mean, Y)

### Plots ----

mydata_c <- data.frame(X_mean = X_mean, Y = Y)

p_c <- ggplot(data = mydata_c, aes(x = X_mean, y = Y)) + 
  geom_point(size=1, alpha=0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Avg. Exemplar Valence", y = "Prosocialness") +
  xlim(c(0,1)) + ylim(c(0,1)) + 
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_c

ggsave(file="Plot_Exemplar-Theory-Simulation.tiff", p_c, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Exemplar-Theory-Simulation.jpeg", p_c, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

Resid2 <- (Y - (0.25+0.5*X_mean))^2

plot(X, Resid2)


mydata_c <- data.frame(Z = Z, Resid2 = Resid2)
p_c_resid2 <- ggplot(data = mydata_c, aes(x = Z, y = Resid2)) + 
  geom_point(size=1, alpha=0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  # geom_density_2d(alpha = 0.5) +
  # geom_coxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.3) +
  labs(x = "Var. Exemplar Valence", y = expression(epsilon^2)) +
  ylim(c(0,0.15)) + xlim(c(0,1)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_c_resid2

ggsave(file="Plot_Exemplar-Theory-Simulation-Residuals.tiff", p_c_resid2, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Exemplar-Theory-Simulation-Residuals.jpeg", p_c_resid2, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

# combined
mydata_c_comb <- data.frame(X_mean = X_mean, Z = Z, Y = Y)
p_c_comb <- ggplot(data = mydata_c_comb, aes(x = X_mean, y = Y, color = Z, fill = Z)) + 
  geom_point(size=1, alpha=0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  # geom_density_2d(alpha = 0.5) +
  # geom_c_comboxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.3) +
  labs(x = "Avg. Exemplar Valence", y = "Prosocialness") +
  ylim(c(0,1)) + 
  scale_color_viridis(begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_c_comb

ggsave(file="Plot_Exemplar-Theory-Simulation_combined.tiff", p_c_comb, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Exemplar-Theory-Simulation_combined.jpeg", p_c_comb, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

### Test Analysis ----

mydata <- data.frame(X = X_mean, Z = Z, Y = Y)

quant <- rq(formula = Y ~ X+Z, tau = c(0.1,0.25,0.5,0.75,0.9), data = mydata)

#### quantile coef plot ----
# OLS
lm <- lm(data=mydata, 
         formula =  Y ~  X+Z)

ols <- as.data.frame(coef(lm))
ols.ci <- as.data.frame(confint(lm))
ols2 <- cbind(ols, ols.ci)
ols2 <- tibble::rownames_to_column(ols2, var="term") |> 
  dplyr::filter(!grepl("(Intercept)", term))


quantile_coef_plot <- quant |> broom::tidy(se.type = "rank", conf.int = TRUE, conf.level = 0.95) |> 
  dplyr::filter(!grepl("(Intercept)", term)) %>%
  ggplot(aes(x=tau,y=estimate))+
  # quantilie results
  geom_point(color="#27408b", size = 1)+ 
  geom_line(color="#27408b")+ 
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.25, fill="#27408b")+
  
  # OLS results
  geom_hline(data = ols2, aes(yintercept= `coef(lm)`), lty=1, color="red")+
  geom_hline(data = ols2, aes(yintercept= `2.5 %`), lty=2, color="red")+
  geom_hline(data = ols2, aes(yintercept= `97.5 %`), lty=2, color="red")+
  scale_x_continuous(breaks=quantiles) +
  ylab("coef") +
  facet_wrap(~term, scales="free", ncol=2)

quantile_coef_plot

ggsave(file="Plot_Prototype-Theory-b-Quantile-coef-plot.tiff", quantile_coef_plot, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Prototype-Theory-b-Quantile-coef-plot.jpeg", quantile_coef_plot, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


summary(quant)
anova(quant, joint = FALSE)


vfr <- dglm(data=mydata,
            formula = Y~X+Z,
            dformula = Y~X+Z, 
            family = gaussian(link = "identity"), dlink = "log", method="reml")
summary(vfr)
anova(vfr)

p_c_quantile <- ggplot(data = mydata, aes(x = Z, y = Y)) + 
  geom_point(alpha=0.3, size=1) +
  #geom_boxplot(aes(group= cut_width(Z,0.05)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.3) +
  geom_quantile(
    method = "rqss", quantiles = c(0.1,0.25,0.5,0.75,0.9)
    #formula = y ~ x, quantiles = c(0.1,0.25,0.5,0.75,0.9)
    , color="darkblue") +
  # geom_density_2d(alpha = 0.5) +
  # geom_boxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.3) +
  labs(x = "Var. Exemplar Valence", y = "Prosocialness") +
  ylim(c(0,1)) + xlim(c(0,1)) +
  scale_color_viridis(begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_c_quantile

ggsave(file="Plot_Exemplar-Theory-Simulation-Quantile-Plot.tiff", p_c_quantile, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Exemplar-Theory-Simulation-Quantile-Plot.jpeg", p_c_quantile, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

out50 <- rq(Y ~ X + Z, 
            data=mydata,tau=c(0.5))
out90 <- rq(Y ~ X + Z, 
            data=mydata,tau=c(0.9))
out75 <- rq(Y ~ X + Z, 
            data=mydata,tau=c(0.75))
out25 <- rq(Y ~ X + Z, 
            data=mydata,tau=c(0.25))
out10 <- rq(Y ~ X + Z, 
            data=mydata,tau=c(0.1))

predictdata <- mydata
predictdata$X = median(mydata$X)


predictdata[,c("yhatq50","yhatq50.lo","yhatq50.hi")] <- predict(out50, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("yhatq90","yhatq90.lo","yhatq90.hi")] <- predict(out90, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("yhatq75","yhatq75.lo","yhatq75.hi")] <- predict(out75, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("yhatq25","yhatq25.lo","yhatq25.hi")] <- predict(out25, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("yhatq10","yhatq10.lo","yhatq10.hi")] <- predict(out10, newdata=predictdata, interval="confidence", level = 0.95)


quant.predict.final <- ggplot(predictdata, aes(x=Z, y=Y)) + 
  geom_line(aes(x=Z, y=yhatq90)) +
  geom_ribbon(aes(ymin=yhatq90.lo, ymax=yhatq90.hi), alpha=0.3, linetype="dotted") +
  geom_line(aes(x=Z, y=yhatq75)) +
  geom_ribbon(aes(ymin=yhatq75.lo, ymax=yhatq75.hi), alpha=0.3, linetype="dotted") +
  geom_line(aes(x=Z, y=yhatq50)) +
  geom_ribbon(aes(ymin=yhatq50.lo, ymax=yhatq50.hi), alpha=0.3, linetype="dotted") +
  geom_line(aes(x=Z, y=yhatq25)) +
  geom_ribbon(aes(ymin=yhatq25.lo, ymax=yhatq25.hi), alpha=0.3, linetype="dotted") +
  geom_line(aes(x=Z, y=yhatq10)) +
  geom_ribbon(aes(ymin=yhatq10.lo, ymax=yhatq10.hi), alpha=0.3, linetype="dotted") +
  ylim(c(0,1)) + xlim(c(0,1)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

quant.predict.final

ggsave(file="Plot_Exemplar-Theory-Simulation-Predicted-Quantile-Plot.tiff", quant.predict.final, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Exemplar-Theory-Simulation-Predicted-Quantile-Plot.jpeg", quant.predict.final, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

# predictdglm <- mydata
# predictdglm$Z <- median(predictdglm$Z)
# predictdglm$yhatv <- predict(vfr, dispersion =1, newdata=predictdglm)
# predictdglm$upperhatv <- predictdglm$yhatv+qnorm(0.75,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
# predictdglm$highesthatv <- predictdglm$yhatv+qnorm(0.9,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
# predictdglm$lowerhatv <- predictdglm$yhatv+qnorm(0.25,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
# predictdglm$lowesthatv <- predictdglm$yhatv+qnorm(0.1,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
# 
# dglm.predict.final <- ggplot(predictdglm, aes(x=X, y=Y)) +
#   geom_line(aes(x=X, y=yhatv)) +
#   geom_line(aes(x=X, y=highesthatv)) +
#   geom_line(aes(x=X, y=upperhatv)) +
#   geom_line(aes(x=X, y=lowerhatv)) +
#   geom_line(aes(x=X, y=lowesthatv)) +
#   xlim(c(0,1)) +
#   scale_color_viridis(begin=0.1, end=0.7, option="inferno") +
#   scale_fill_viridis(begin=0.1, end=0.7, option="inferno") +
#   guides(shape = guide_legend(override.aes = list(alpha= 1)))
# 
# dglm.predict.final

predictdglm <- mydata
predictdglm$X <- median(predictdglm$X)
predictdglm$yhatv <- predict(vfr, dispersion =1, newdata=predictdglm)
predictdglm$upperhatv <- predictdglm$yhatv+qnorm(0.75,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
predictdglm$highesthatv <- predictdglm$yhatv+qnorm(0.9,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
predictdglm$lowerhatv <- predictdglm$yhatv+qnorm(0.25,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
predictdglm$lowesthatv <- predictdglm$yhatv+qnorm(0.1,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))


logvarpredict <- predict(vfr$dispersion.fit, dispersion=2, newdata=predictdglm, se.fit=TRUE)


predictdglm[,"sdhat"] <- sqrt(exp(logvarpredict$fit))
predictdglm[,"sdlo"] <- sqrt(exp(logvarpredict$fit - qnorm(p=0.975,0,1)*logvarpredict$se.fit))
predictdglm[,"sdhi"] <- sqrt(exp(logvarpredict$fit + qnorm(p=0.975,0,1)*logvarpredict$se.fit))

sd.predict.final <- ggplot(predictdglm, aes(x=Z, y=sdhat)) +
  geom_line(aes(x=Z, y=sdhat)) +
  geom_ribbon(aes(ymin=sdlo, ymax=sdhi), alpha=0.3, linetype="dotted") +
  labs(x="Var. Exemplar Valence",y=bquote(widehat(sigma)(Prosocialness)))+
  xlim(c(0,1)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

sd.predict.final

ggsave(file="Plot_Exemplar-Theory-Simulation-Predicted-SD-Plot.tiff", sd.predict.final, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Exemplar-Theory-b-Simulation-Predicted-SD-Plot.jpeg", sd.predict.final, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

## Combine plots ----

p_b_both <- p_b + p_b_resid2 +
  plot_layout(nrow=1, ncol = 2)

p_b_both

ggsave(file="Plot_Simulation_Data_1b_together.tiff", p_b_both, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Simulation_Data_1b_together.jpeg", p_b_both, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

p_c_both <- p_c + p_c_resid2 +
  plot_layout(nrow=1, ncol = 2)

p_c_both

ggsave(file="Plot_Simulation_Data_2_together.tiff", p_c_both, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Simulation_Data_2_together.jpeg", p_c_both, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

p_b_c <- ((p_b | p_c) + plot_layout(axes = "collect")) /
  ((p_b_resid2 | p_c_resid2) +   plot_layout(axes = "collect"))

p_b_c <- p_b_c + plot_annotation(tag_levels = "a",
                                       tag_suffix = ")",
                                       tag_prefix = "(") &
theme(plot.tag.position = c("topleft"),
      plot.tag.location = c("panel"),
        plot.tag = element_text(size = 8, hjust = 0.5, vjust = 0, color = "black", face = "bold", family = "sans")) 

p_b_c

ggsave(file="Plot_Simulation_Data_1b_2_together.tiff", p_b_c, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Simulation_Data_1b_2_together.jpeg", p_b_c, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


p_all <- ((p_a | p_b | p_c) + plot_layout(axes = "collect_y")) /
  ((p_a_resid2 | p_b_resid2 | p_c_resid2) +   plot_layout(axes = "collect_y"))

p_all <- p_all + plot_annotation(tag_levels = "a",
                                 tag_suffix = ")",
                                 tag_prefix = "(") &
  theme(plot.tag.position = c("topleft"),
        plot.tag.location = c("panel"),
        plot.tag = element_text(size = 8, hjust = 0.5, vjust = 0, color = "black", face = "bold", family = "sans")) 

p_all

ggsave(file="Plot_Simulation_Data_all_together.tiff", p_all, 
       device = tiff, width = 114.3, height = 1.5*88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="Plot_Simulation_Data_all_together.jpeg", p_all, 
       device = jpeg, width = 114.3, height = 1.5*88.9, units = "mm", limitsize = FALSE, dpi = 600)