##############HEADER##############################
## author:    Daniel Gotthardt
## contact:   daniel.gotthardt@uni-hamburg.de / daniel.gotthardt@gmx.de
## file name: Alice_Empirical_Example.R
## Context:   Theorizing and Testing Variability in Social Cognition
## Input:     -
## Output:    Variance and Quantile coefficient and prediction plots
## Summary:   This R-File analyses the generated empirical data
# with quantile regression and variance function regression

################BODY##############################

# Setup ----------------------
## Clean up workspace -------------------------------------
rm(list=ls(all.names = TRUE))
gc()

## Load packages -----
library(ggplot2)
library(viridis)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(dglm)
library(quantreg)
library(broom)
library(texreg)
library(marginaleffects)

## Load data -----

mydata <- readRDS("data/Simulated_Empirical_Data.Rds")

# Graphical exploration ----

theme_set(theme_bw(base_size = 10,
                   base_family = "sans"))


p <- ggplot(data = mydata, aes(x = X.mean, y = Y.obs)) + 
  geom_point(size=2, alpha=0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Valence", y = "Prosocialness") +
  ylim(c(0,1)) + xlim(c(0,1)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p

# OLS
lm <- lm(data=mydata, 
         formula =  Y.obs ~  X.mean+I(X.mean^2)+Z)

mydata$Y_pred <- predict(lm)

mydata$Resid <- sqrt((mydata$Y.obs - mydata$Y_pred)^2)

p_b_resid <- ggplot(data = mydata, aes(x = X.mean, y = Resid)) + 
  geom_point(size=2, alpha=0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  # geom_density_2d(alpha = 0.5) +
  # geom_boxplot(aes(group= cut_width(X,0.2)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.5) +
  ylim(c(0,0.4)) + xlim(c(0,1)) +
  labs(x = "Valence", y = expression(sqrt(epsilon^2))) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

p_b_resid

# Model fitting

## Fit quantreg ----

quantiles <- c(0.1,0.25,0.5,0.75,0.9)

quant <- rq(formula = Y.obs ~ X.mean+I(X.mean^2)+Z, 
            tau = quantiles, 
            data = mydata)
summary(quant)

quant_ano_full <- anova(quant, joint = FALSE, iid = FALSE)


out50 <- rq(Y.obs ~ X.mean+I(X.mean^2)+Z, 
            data=mydata,tau=c(0.5))
out90 <- rq(Y.obs ~ X.mean+I(X.mean^2)+Z, 
            data=mydata,tau=c(0.9))
out75 <- rq(Y.obs ~ X.mean+I(X.mean^2)+Z, 
            data=mydata,tau=c(0.75))
out25 <- rq(Y.obs ~ X.mean+I(X.mean^2)+Z, 
            data=mydata,tau=c(0.25))
out10 <- rq(Y.obs ~ X.mean+I(X.mean^2)+Z, 
            data=mydata,tau=c(0.1))


quantlist = list(tau10 = out10, tau25 = out25, tau50 = out50, tau75 = out75, tau90 = out90)

quantanova <- lapply(quantlist, function(m1) lapply(quantlist[!(quantlist %in% list(m1))], function(m2) {
  anova(m1, m2, iid= FALSE, joint = FALSE)
}
)
)

ano2575 <- quantanova$tau25$tau75
ano1090 <- quantanova$tau10$tau90
ano1090

screenreg(quantlist, digits=2,
split.table = Inf)





### quantile coef plot ----

coef_labels <- list(
  'I(X.mean^2)'=expression(X^2),
  'X.mean'="X"
  # ,
  # "Z" = "Z"
)

coef_labeller <- function(variable,value){
  return(coef_labels[value])
}

ols <- as.data.frame(coef(lm))
ols.ci <- as.data.frame(confint(lm))
ols2 <- cbind(ols, ols.ci)
ols2 <- tibble::rownames_to_column(ols2, var="term") |> 
  dplyr::filter(!grepl("(Intercept)|Z", term))

quantile_coef_plot <- quant |> broom::tidy(se.type = "rank", conf.int = TRUE, conf.level = 0.95) |> 
  dplyr::filter(!grepl("(Intercept)|Z", term)) %>%
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

ggsave(file="plots/Plot_Empirical-Quantile-Coef-Plot.tiff", quantile_coef_plot,
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/Plot_Empirical-Quantile-Coef-Plot.jpeg", quantile_coef_plot,
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

summary(quant)
anova(quant, joint = FALSE)

## Fit dglm ----

vfr <- dglm(data=mydata,
            formula = Y.obs ~ X.mean+I(X.mean^2)+Z,
            dformula = Y.obs ~ X.mean+I(X.mean^2)+Z, 
            family = gaussian(link = "identity"), 
            dlink = "log", 
            method="reml")
summary(vfr)

slopes(vfr, variables = "X.mean", dispersion =1)
plot_slopes(vfr, variables = "X.mean",
            condition = "X.mean", dispersion =1)
avg_slopes(vfr, variables = "X.mean", dispersion =1)


plot_slopes(vfr$dispersion.fit, variables = "X.mean",
            condition = "X.mean", dispersion =2)
avg_slopes(vfr$dispersion.fit, variables = "X.mean", dispersion =2)



vfr_table <- list(beta = vfr, 
                  lambda = vfr$dispersion.fit)

screenreg(vfr_table, dispersion=c(1,2),
          include.aic=FALSE,
          include.bic=FALSE, 
          digits=2
        #file = "vfr.doc"
        )


# Plot predictions ----

## Predictions from qr ----




predictdata <- mydata
predictdata$Z <- median(mydata$Z)
predictdata[,c("q50_yhat","q50_lo","q50_hi")] <- predict(out50, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("q90_yhat","q90_lo","q90_hi")] <- predict(out90, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("q75_yhat","q75_lo","q75_hi")] <- predict(out75, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("q25_yhat","q25_lo","q25_hi")] <- predict(out25, newdata=predictdata, interval="confidence", level = 0.95)
predictdata[,c("q10_yhat","q10_lo","q10_hi")] <- predict(out10, newdata=predictdata, interval="confidence", level = 0.95)

quant.predict.final <- predictdata |> 
  tidyr::pivot_longer(starts_with("q"),
                      names_to = c("tau", ".value"),
                      names_pattern = "^(.*)_(.*)$") |> 
  ggplot(aes(x=X.mean, y=Y.obs, color = tau, fill = tau, group = tau)) + 
  geom_line(aes(x=X.mean, y=yhat)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.3, linetype="dotted") +
  ylim(c(-0.25,1.00)) + xlim(c(0,1)) +
  labs(x="Stereotype Valence", y = bquote(widehat(Q)[tau](Prosocialness))) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno",
                      name = "Quantiles",
                      labels = quantiles) +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno",
                     name = "Quantiles",
                     labels = quantiles) +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

quant.predict.final

## Predictions from dglm ----

ggsave(file="plots/Plot_Empirical-Predicted-Quantile-Plot.tiff", quant.predict.final,
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/Plot_Empirical-Predicted-Quantile-Plot.jpeg", quant.predict.final,
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

predictdglm <- mydata
predictdglm$Z <- median(mydata$Z)

meanpredict <- predict(vfr, dispersion = 1, newdata = predictdglm, se.fit=TRUE)

predictdglm[,"yhat"] <- meanpredict$fit
predictdglm[,"ci_lo"] <- meanpredict$fit - qnorm(p=0.975,0,1)*meanpredict$se.fit
predictdglm[,"ci_hi"] <- meanpredict$fit + qnorm(p=0.975,0,1)*meanpredict$se.fit

predictdglm$predquart_hi <- predictdglm$yhat+qnorm(0.75,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
predictdglm$predperc_hi <- predictdglm$yhat+qnorm(0.9,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
predictdglm$predquart_lo <- predictdglm$yhat+qnorm(0.25,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))
predictdglm$predperc_lo <- predictdglm$yhat+qnorm(0.1,0,1)*sqrt(predict(vfr$dispersion.fit, dispersion=2, type ="response", newdata=predictdglm))

dglm.predict.final <- predictdglm  |>  
#   tidyr::pivot_longer(contains("_"),
#                       names_to = c("interval", ".value"),
#                       names_pattern = "^(.*)_(.*)$") |> 
#   mutate(interval = factor(interval, levels = c("predperc", "predquart","ci"))) |> 
  ggplot(aes(x=X.mean, y=Y.obs)) +
  geom_line(aes(x=X.mean, y=yhat)) +
  geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi, 
                  #group = interval, fill = interval, color = interval, linetype = interval
                  ), 
              alpha=0.4) +
  labs(y=bquote(widehat(E)(Prosocialness)), x = "Stereotype Valence")+
  xlim(c(0,1)) + ylim(c(0,1)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

dglm.predict.final

ggsave(file="plots/Plot_Empirical-Predicted-Mean-Plot.tiff", dglm.predict.final,
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/Plot_Empirical-Predicted-Mean-Plot.jpeg", dglm.predict.final,
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

logvarpredict <- predict(vfr$dispersion.fit, dispersion=2, newdata=predictdglm, se.fit=TRUE)

predictdglm[,"sdhat"] <- sqrt(exp(logvarpredict$fit))
predictdglm[,"sdlo"] <- sqrt(exp(logvarpredict$fit - qnorm(p=0.975,0,1)*logvarpredict$se.fit))
predictdglm[,"sdhi"] <- sqrt(exp(logvarpredict$fit + qnorm(p=0.975,0,1)*logvarpredict$se.fit))

sd.predict.final <- ggplot(predictdglm, aes(x=X.mean, y=sdhat)) +
  geom_line(aes(x=X.mean, y=sdhat)) +
  geom_ribbon(aes(ymin=sdlo, ymax=sdhi), alpha=0.3, linetype="dotted") +
  labs(y=bquote(widehat(sigma)(Prosocialness)), x = "Stereotype Valence")+
  xlim(c(0,1)) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))

sd.predict.final

ggsave(file="plots/Plot_Empirical-Predicted-SD-Plot.tiff", sd.predict.final,
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/Plot_Empirical-Predicted-SD-Plot.jpeg", sd.predict.final,
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

p_vfr <- dglm.predict.final + sd.predict.final +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")",
                  tag_prefix = "(") +
  plot_layout(nrow=1, ncol = 2,
              axis_titles = "collect_x") &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 8, hjust = -1, vjust = 0, color = "black", face = "bold", family = "sans")) 

p_vfr

ggsave(file="plots/Plot_Empirical-Predicted-VFR-plot.tiff", p_vfr, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/Plot_Empirical-Predicted-VFR-plot.jpeg", p_vfr, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)