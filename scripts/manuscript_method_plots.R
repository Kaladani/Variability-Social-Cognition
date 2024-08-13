

# Clean up workspace -------------------------------------
rm(list=ls(all.names = TRUE))
gc()


#install.packages("ggplot2")
#install.packages("viridis")
library(ggplot2)
library(viridis)
library(tidyverse)

theme_set(theme_bw(base_size = 10,
                   base_family = "sans"))


set.seed(4242)

X <- runif(100, 1, 10)
normal <- rnorm(100,mean=0,sd=1)
u <- 0.5*X*normal
Y <- X + u
dat <- data.frame(X,Y)

lm.fit <- lm(Y~X, data=dat)
dat$resid <- residuals(lm.fit)
dat$resid2 <- dat$resid^2


# plots for conditional distribution

scatter <- ggplot(dat, aes(x=X, y=Y)) + 
  geom_point(alpha=0.5) +
  geom_smooth(formula=y~x, method="lm", linetype="dashed", color="red", fill=NA)


ggsave(file="plots/scatter.tiff", scatter, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/scatter.jpeg", scatter, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


resid.sq <- ggplot(dat, aes(x=X, y=resid2)) + 
  geom_point(alpha=0.5) + 
  geom_hline(aes(yintercept=mean(resid2)), linetype="dashed", color="red", size=1) +
  stat_summary_bin(fun="mean",
                   binwidth=1,
                   color="blue",
                   geom="segment",
                   aes(xend=..x..-0.5,
                       yend=..y..),
                   size=1) +
  stat_summary_bin(fun="mean",
                   binwidth=1,
                   color="blue",
                   geom="segment",
                   aes(xend=..x..+0.5,
                       yend=..y..),
                   size=1) +
  labs(y=bquote(epsilon^2))

ggsave(file="plots/resid.sq.tiff", resid.sq, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/resid.sq.jpeg", resid.sq, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


scatterbox <- ggplot(dat, aes(x=X, y=Y)) + 
  geom_point(alpha=0.5) +
  geom_boxplot(aes(group= cut_width(X,1)), outlier.shape=NA, varwidth=TRUE,fill="grey", color="black", alpha=0.5) +
  geom_quantile(formula=y~1,quantiles=c(0.25,0.75), size=1, linetype="dashed", color="red")   

ggsave(file="plots/scatterbox.tiff", scatterbox, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/scatterbox.jpeg", scatterbox, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


## Plots for methods section

scatterbox_quant <- ggplot(dat, aes(x=X, y=Y)) + 
  geom_point(alpha=0.5, size=1) +
  geom_boxplot(aes(group= cut_width(X,1)), outlier.shape=NA, varwidth=TRUE,fill="grey", alpha=0.5) +
  geom_quantile(formula=y~x,quantiles=c(0.25,0.75), size=1, color="darkblue") 

ggsave(file="plots/scatterbox_quant.tiff", scatterbox_quant, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/scatterbox_quant.jpeg", scatterbox_quant, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

resid.sq_pred <- ggplot(dat, aes(x=X, y=resid2)) + 
  geom_point(alpha=0.5)  +
  geom_smooth(formula=y~x, method="glm", method.args = list(family= Gamma (link=log)), fill=NA, color="darkblue",
              size=1) +
  geom_smooth(formula=y~x, method="lm", linetype="longdash", fill=NA, color="orange") +
  labs(y=bquote(epsilon^2))

ggsave(file="plots/resid.sq_pred.tiff", resid.sq_pred, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/resid.sq_pred.jpeg", resid.sq_pred, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


#Define convenience functions

## Check for extreme outliers
is.extreme.Fun <- function(x) {
  return(x < quantile(x, 0.25, na.rm=TRUE) - 3 * IQR(x, na.rm=TRUE) | x > quantile(x, 0.75, na.rm=TRUE) + 3 * IQR(x, na.rm=TRUE))
}


set.seed(42)

## Plots for Simulation

## Simulate exponential variance

n_obs.1000 <- 1000

# Set nonrandom parameters

beta_0 <- 3
beta_1 <- 0.5
beta_2 <- 0.2
beta_3 <- -1
beta_4 <- -3

lambda_0 <- 2
lambda_1 <- 0.5
lambda_2 <- 0.2
lambda_3 <- -0.4
sigma <- 1

# Draw random data

x <- runif(n_obs.1000, 1, 10)
xx <- x^2
z <- rbinom(n_obs.1000,1,0.5)
xz <- x*z
u <- sqrt(exp(lambda_0 + lambda_1*x + lambda_2*z + lambda_3*x*z))*rnorm(n_obs.1000,mean=0,sd=sigma)
y <- beta_0+x*beta_1 + x^2*beta_2 + z*beta_3 + x*z*beta_4 + u
y.true <- beta_0+x*beta_1 + x^2*beta_2 + z*beta_3 + x*z*beta_4

expdf <- data.frame(x,xx,z,xz,u,y,y.true)


## Simulate quadratic variance

n_obs.1000 <- 1000

# Set nonrandom parameters

beta_0 <- 3
beta_1 <- 0.5
beta_2 <- 0.2
beta_3 <- -1
beta_4 <- -3

lambda_0 <- 5
lambda_1 <- 3
lambda_2 <- 1
lambda_3 <- -2
sigma <- 1

# Draw random data
x <- runif(n_obs.1000, 1, 10)
xx <- x^2
z <- rbinom(n_obs.1000,1,0.5)
xz <- x*z
u <- (lambda_0+lambda_1*x+lambda_2*z+lambda_3*x*z)*rnorm(n_obs.1000,mean=0,sd=sigma)
y <- beta_0+x*beta_1 + x^2*beta_2 + z*beta_3 + x*z*beta_4 + u
y.true <- beta_0+x*beta_1 + x^2*beta_2 + z*beta_3 + x*z*beta_4

lindf <- data.frame(x,xx,z,xz,u,y,y.true)



# Create plot

simdf <- rbind(data.frame(expdf, Scenario=rep("Scenario 6",1000)),
               data.frame(lindf, Scenario=rep("Scenario 7",1000)))
simdf$z <- factor(simdf$z, labels=c("z=0","z=1"))


simexample <- ggplot(simdf, aes(x=x,y=y, color=Scenario, linetype=z)) + 
  geom_point(size=0.2, alpha=0.5) +
  geom_line(aes(x=x, y=y.true), size=1) +
  facet_grid(z~Scenario) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  scale_fill_viridis(discrete = TRUE, begin=0.1, end=0.7, option="inferno") +
  guides(shape = guide_legend(override.aes = list(alpha= 1)))




simexample

ggsave(file="plots/simexample.tiff", simexample, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/simexample.jpeg", simexample, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


# I Scenario Null with Non-linear mean model --------------------------
## Load Data

est.long.30 <- readRDs(file="data/S3_Nullx2_est_long_30.rds")
est.30.agg <- readRDs(file="data/S3_Nullx2_est_agg_30.rds")
est.long.50 <- readRDs(file="data/S3_Nullx2_est_long_50.rds")
est.50.agg <- readRDs(file="data/S3_Nullx2_est_agg_50.rds")
est.long.1000 <- readRDs(file="data/S3_Nullx2_est_long_1000.rds")
est.1000.agg <- readRDs(file="data/S3_Nullx2_est_agg_1000.rds")

## Data preperation

est.n <- merge(est.long.30,est.long.50, by=c("repetition","model","coef.type"), suffixes=c(".30",".50"))
est.n <- merge(est.n, est.long.1000, by=c("repetition","model","coef.type"))
colnames(est.n)[(ncol(est.n)-8):ncol(est.n)] <- c( "coef.1000", "se.1000", "p.1000","true.1000", "bias.1000", "ci_lo.1000", "ci_hi.1000", "cover.1000", "reject.1000")
est.n.graph <- reshape(est.n, 
                       direction ="long", 
                       idvar= c("repetition","model","coef.type"), 
                       varying = 4:(ncol(est.n)),
                       timevar ="n_obs")
est.n.graph$coef.type <- sapply(strsplit(est.n.graph$coef.type,"_"),
                                function(vec) 
                                  do.call(sprintf, c(list("widehat(%s)[%s]"),vec
                                  )))
est.n.graph$n_obs <- factor(sprintf("n = %s",est.n.graph$n_obs),
                            levels=c("n = 30", "n = 50", "n = 1000"))
est.n.graph$model <- factor(est.n.graph$model,
                            levels=c("Two-Step", "ML", "REML","QR"))

est.n.agg <- merge(est.30.agg,est.50.agg, by=c("repetition","model","coef.type"), suffixes=c(".30",".50"))
est.n.agg <- merge(est.n.agg, est.1000.agg, by=c("repetition","model","coef.type"))
colnames(est.n.agg)[(ncol(est.n.agg)-12):ncol(est.n.agg)] <- sprintf("%s.1000",colnames(est.n.agg)[(ncol(est.n.agg)-12):ncol(est.n.agg)])
est.n.agg.graph <- reshape(est.n.agg, 
                           direction ="long", 
                           idvar= c("repetition","model","coef.type"), 
                           varying = 4:(ncol(est.n.agg)),
                           timevar ="n_obs")
est.n.agg.graph$coef.type <- sapply(strsplit(est.n.agg.graph$coef.type,"_"),
                                    function(vec) 
                                      do.call(sprintf, c(list("widehat(%s)[%s]"),vec
                                      )))
est.n.agg.graph$n_obs <- factor(sprintf("n = %s",est.n.agg.graph$n_obs),
                                levels=c("n = 30", "n = 50", "n = 1000"))
est.n.agg.graph$model <- factor(est.n.agg.graph$model,
                                levels=c("Two-Step", "ML", "REML","QR"))

## Delete temporary variables
rm(list=ls()[! (ls() %in% lsf.str() |ls() %in% c("est.n.agg.graph","est.n.graph"))])



## Histogram of lambda_1 for all n_obs ####

hist.coef.n <- ggplot(data=est.n.graph, aes(x=coef)) +
  geom_histogram(color="black", fill="grey",binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T) + #Freedman-Diaconis Rule for binwidth
  labs(x="",y="f")

hist.coef.n %+% subset(est.n.graph, coef.type %in% "widehat(lambda)[1]") +
  facet_grid(n_obs~model) + 
  geom_vline(data=subset(est.n.agg.graph, coef.type %in% "widehat(lambda)[1]"), aes(xintercept=coef),
             size=1, linetype="dashed") +
  labs(x=bquote(widehat(lambda)[1]))                                           


## Boxplot of lambda_1 for all n_obs ####

box.coef.n <- ggplot(data=est.n.graph, aes(y=model, x=coef)) +
  geom_boxplot(width=0.5)  

box.coef.split <- box.coef.n %+% subset(est.n.graph, coef.type %in% "widehat(lambda)[1]") +
  facet_grid(n_obs~.) + 
  labs(x=bquote(widehat(lambda)[1]), y="method") + 
  geom_vline(data=subset(est.n.agg.graph, coef.type %in% "widehat(lambda)[1]"), aes(xintercept=0),
             size=0.5, color="grey", linetype="dashed")

ggsave(file="plots/S1-misspecified-boxplot.tiff", box.coef.split, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/S1-misspecified-boxplot.jpeg", box.coef.split, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


## Check for outliers

est.n.graph.split <- split(est.n.graph, list(est.n.graph$n_obs, est.n.graph$coef.type, est.n.graph$model))
est.n.graph.split <- lapply(est.n.graph.split, transform, outlier = is.extreme.Fun(coef))
est.n.graph <- unsplit(est.n.graph.split,list(est.n.graph$n_obs, est.n.graph$coef.type, est.n.graph$model))

est.n.graph[which(est.n.graph$outlier),c("n_obs","coef.type","model","repetition","coef")]
est.n.graph[which(est.n.graph$outlier),c("coef")] <- NA

## Histogram without outlier

hist.coef.n <- ggplot(data=est.n.graph, aes(x=coef)) +
  geom_histogram(color="black", fill="grey",binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T) + #Freedman-Diaconis Rule for binwidth
  labs(x="",y="f")

hist.coef.split <- hist.coef.n %+% subset(est.n.graph, coef.type %in% "widehat(lambda)[1]") +
  facet_grid(n_obs~model) +
  geom_vline(data=subset(est.n.agg.graph, coef.type %in% "widehat(lambda)[1]"), aes(xintercept=coef),
             size=1, linetype="dashed") +
  labs(x=bquote(widehat(lambda)[1]))


ggsave(file="plots/S1-misspecified-histogram.tiff", hist.coef.split, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/S1-misspecified-histogram.jpeg", hist.coef.split, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


# II Scenario Null with vertical Outlier --------------------------
## Load Data

est.long.30 <- readRDs(file="data/S4_Nullout_est_long_30.rds")
est.30.agg <- readRDs(file="data/S4_Nullout_est_agg_30.rds")
est.long.50 <- readRDs(file="data/S4_Nullout_est_long_50.rds")
est.50.agg <- readRDs(file="data/S4_Nullout_est_agg_50.rds")
est.long.1000 <- readRDs(file="data/S4_Nullout_est_long_1000.rds")
est.1000.agg <- readRDs(file="data/S4_Nullout_est_agg_1000.rds")

## Data preperation

est.n <- merge(est.long.30,est.long.50, by=c("repetition","model","coef.type"), suffixes=c(".30",".50"))
est.n <- merge(est.n, est.long.1000, by=c("repetition","model","coef.type"))
colnames(est.n)[(ncol(est.n)-8):ncol(est.n)] <- c( "coef.1000", "se.1000", "p.1000","true.1000", "bias.1000", "ci_lo.1000", "ci_hi.1000", "cover.1000", "reject.1000")
est.n.graph <- reshape(est.n, 
                       direction ="long", 
                       idvar= c("repetition","model","coef.type"), 
                       varying = 4:(ncol(est.n)),
                       timevar ="n_obs")
est.n.graph$coef.type <- sapply(strsplit(est.n.graph$coef.type,"_"),
                                function(vec) 
                                  do.call(sprintf, c(list("widehat(%s)[%s]"),vec
                                  )))
est.n.graph$n_obs <- factor(sprintf("n = %s",est.n.graph$n_obs),
                            levels=c("n = 30", "n = 50", "n = 1000"))
est.n.graph$model <- factor(est.n.graph$model,
                            levels=c("Two-Step", "ML", "REML","QR"))

est.n.agg <- merge(est.30.agg,est.50.agg, by=c("repetition","model","coef.type"), suffixes=c(".30",".50"))
est.n.agg <- merge(est.n.agg, est.1000.agg, by=c("repetition","model","coef.type"))
colnames(est.n.agg)[(ncol(est.n.agg)-12):ncol(est.n.agg)] <- sprintf("%s.1000",colnames(est.n.agg)[(ncol(est.n.agg)-12):ncol(est.n.agg)])
est.n.agg.graph <- reshape(est.n.agg, 
                           direction ="long", 
                           idvar= c("repetition","model","coef.type"), 
                           varying = 4:(ncol(est.n.agg)),
                           timevar ="n_obs")
est.n.agg.graph$coef.type <- sapply(strsplit(est.n.agg.graph$coef.type,"_"),
                                    function(vec) 
                                      do.call(sprintf, c(list("widehat(%s)[%s]"),vec
                                      )))
est.n.agg.graph$n_obs <- factor(sprintf("n = %s",est.n.agg.graph$n_obs),
                                levels=c("n = 30", "n = 50", "n = 1000"))
est.n.agg.graph$model <- factor(est.n.agg.graph$model,
                                levels=c("Two-Step", "ML", "REML","QR"))

## Delete temporary variables
rm(list=ls()[! (ls() %in% lsf.str() |ls() %in% c("est.n.agg.graph","est.n.graph"))])

## Check for outliers

est.n.graph.split <- split(est.n.graph, list(est.n.graph$n_obs, est.n.graph$coef.type, est.n.graph$model))
est.n.graph.split <- lapply(est.n.graph.split, transform, outlier = is.extreme.Fun(coef))
est.n.graph <- unsplit(est.n.graph.split,list(est.n.graph$n_obs, est.n.graph$coef.type, est.n.graph$model))

est.n.graph[which(est.n.graph$outlier),c("n_obs","coef.type","model","repetition","coef","se")]

# the ML-fit for repetition 328 is extremely bad (-10^14 bzw. 10^15 statt 0) 
# and needs to be eliminated to allow for plotting; mean is also highly distorted

est.n.graph[which(abs(est.n.graph$coef) >= 100),c("coef","se")] <- NA
est.n.agg.graph[which(abs(est.n.agg.graph$coef) >= 100),c("coef","se")] <- NA


## Histogram of lambda_1 for all n_obs ####

hist.coef.n <- ggplot(data=est.n.graph, aes(x=coef)) +
  geom_histogram(color="black", fill="grey",binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T) + #Freedman-Diaconis Rule for binwidth
  labs(x="",y="f")

hist.coef.split <- hist.coef.n %+% subset(est.n.graph, coef.type %in% "widehat(lambda)[1]") +
  facet_grid(n_obs~model) +
  geom_vline(data=subset(est.n.agg.graph, coef.type %in% "widehat(lambda)[1]"), aes(xintercept=coef),
             size=1, linetype="dashed") +
  labs(x=bquote(widehat(lambda)[1]))

ggsave(file="plots/verticaloutliers-histogram.tiff", hist.coef.split, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/verticaloutliers-histogram.jpeg", hist.coef.split, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


## Boxplot of lambda_1 for all n_obs ####

box.coef.n <- ggplot(data=est.n.graph, aes(y=model, x=coef)) +
  geom_boxplot(width=0.5)  

box.coef.split <- box.coef.n %+% subset(est.n.graph, coef.type %in% "widehat(lambda)[1]") +
  facet_grid(n_obs~.) + 
  labs(x=bquote(widehat(lambda)[1]), y="method") + 
  geom_vline(data=subset(est.n.agg.graph, coef.type %in% "widehat(lambda)[1]"), aes(xintercept=0),
             size=0.5, color="grey", linetype="dashed")

ggsave(file="plots/verticaloutliers-boxplot.tiff", box.coef.split, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/verticaloutliers-boxplot.jpeg", box.coef.split, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)


# ## Histogram of se lambda_1 for all n_obs ####
# 
# hist.se.n <- ggplot(data=est.n.graph, aes(x=se)) +
#   geom_histogram(color="black", fill="grey",binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T) + #Freedman-Diaconis Rule for binwidth
#   labs(x="",y="f")
# 
# hist.se.n %+% subset(est.n.graph, coef.type %in% "widehat(lambda)[1]" & model != "QR") +
#   facet_grid(n_obs~model) + 
#   geom_vline(data=subset(est.n.agg.graph, coef.type %in% "widehat(lambda)[1]" & model != "QR"), aes(xintercept=se),
#              size=1, linetype="dashed") +
#   labs(x=bquote(SE(widehat(lambda)[1])))                                            
# 
# 
# ## Boxplot of se lambda_1 for all n_obs ####
# 
# box.se.n <- ggplot(data=est.n.graph, aes(y=model, x=se)) +
#   geom_boxplot(width=0.5)  
# 
# box.se.n %+% subset(est.n.graph, (coef.type %in% "widehat(lambda)[1]" & model != "QR")) +
#   facet_grid(n_obs~.) + 
#   labs(x=bquote(SE(widehat(lambda)[1])), y="method") + 
#   geom_vline(data=subset(est.n.agg.graph, coef.type %in% "widehat(lambda)[1]" & model != "QR"), aes(xintercept=0),
#              size=0.5, color="grey", linetype="dashed")


## Relative %-Error in ModSE

est.n.agg.graph$relse <- 100*(est.n.agg.graph$modelse/est.n.agg.graph$empse-1)



# III Scenario 7 Lin ----------------------------------------

## Load Data

est.long.30 <- readRDs(file="data/S7_Lin_est_long_30.rds")
est.30.agg <- readRDs(file="data/S7_Lin_est_agg_long_30.rds")
est.long.50 <- readRDs(file="data/S7_Lin_est_long_50.rds")
est.50.agg <- readRDs(file="data/S7_Lin_est_agg_long_50.rds")
est.long.1000 <- readRDs(file="data/S7_Lin_est_long_1000.rds")
est.1000.agg <- readRDs(file="data/S7_Lin_est_agg_long_1000.rds")

## Prepare Data

## Data preperation

est.n <- merge(est.long.30,est.long.50, by=c("repetition","tau","coef.type"), suffixes=c(".30",".50"))
est.n <- merge(est.n, est.long.1000, by=c("repetition","tau","coef.type"))
colnames(est.n)[(ncol(est.n)-4):ncol(est.n)] <- c( "coef.1000", "se.1000", "true.1000", "ci_lo.1000", "ci_hi.1000")
est.n.graph <- reshape(est.n, 
                       direction ="long", 
                       idvar= c("repetition","tau","coef.type"), 
                       varying = 4:(ncol(est.n)),
                       timevar ="n_obs")
# est.n.graph$coef.type <- sapply(strsplit(est.n.graph$coef.type,"_"),
#                                 function(vec) 
#                                   do.call(sprintf, c(list("widehat(%s)[%s]"),vec
#                                   )))
est.n.graph$n_obs <- factor(sprintf("n = %s",est.n.graph$n_obs),
                            levels=c("n = 30", "n = 50", "n = 1000"))
# est.n.graph$model <- factor(est.n.graph$model,
#                             levels=c("Two-Step", "ML", "REML","QR"))


## Delete temporary variables
rm(list=ls()[! (ls() %in% lsf.str() |ls() %in% c("est.n.agg.graph","est.n.graph"))])

## Check for outliers

est.n.graph.split <- split(est.n.graph, list(est.n.graph$n_obs, est.n.graph$coef.type, est.n.graph$tau))
est.n.graph.split <- lapply(est.n.graph.split, transform, outlier = is.extreme.Fun(coef))
est.n.graph <- unsplit(est.n.graph.split,list(est.n.graph$n_obs, est.n.graph$coef.type, est.n.graph$tau))

est.n.graph[which(est.n.graph$outlier),c("n_obs","coef.type","tau","repetition","coef","se")]


# ## Histogram of beta_1 for all n_obs ####
# 
# hist.coef.n <- ggplot(data=est.n.graph, aes(x=coef)) +
#   geom_histogram(color="black", fill="grey",binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T) + #Freedman-Diaconis Rule for binwidth
#   labs(x="",y="f")
# 
# hist.coef.n %+% subset(est.n.graph, coef.type %in% "x") +
#   facet_grid(n_obs~tau) + 
#   geom_vline(aes(xintercept=true),
#               size=1, color="blue", linetype="solid") +
#   labs(x=bquote(widehat(beta)[1]))
# 
# ## Boxplot of beta_1 for all n_obs ####
# 
# box.coef.n <- ggplot(data=est.n.graph, aes(y=tau, x=coef)) +
#   geom_boxplot(width=0.5)  
# 
# box.coef.n %+% subset(est.n.graph, coef.type %in% "x") +
#   facet_grid(n_obs~.) + 
#   labs(x=bquote(widehat(beta)[1]), y="Quantil") + 
#   geom_vline(aes(xintercept=true),
#             size=0.5, color="blue", linetype="solid")

# Calculate diff for each simulation
est.n.diff <- reshape(est.n.graph, direction = "wide",
                      idvar=c("repetition","coef.type","n_obs"),
                      timevar="tau")
est.n.diff$diff <- est.n.diff$coef.Tau75-est.n.diff$coef.Tau25
est.n.diff$truediff <- est.n.diff$true.Tau75-est.n.diff$true.Tau25


## Histogram of beta_1 diff ####

hist.coef.n <- ggplot(data=est.n.diff, aes(x=diff, color=n_obs, fill=n_obs)) +
  geom_histogram(alpha=0.4, position="identity",
                 binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)),na.rm=T) + #Freedman-Diaconis Rule for binwidth
  labs(x="",y="f")  +
  scale_color_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1) +
  scale_fill_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1)

hist.coef.split <- hist.coef.n %+% subset(est.n.diff, (coef.type %in% "x")) +
  geom_vline(aes(xintercept=truediff),
             size=1, color="black", linetype="dashed") +
  labs(x=bquote(widehat(beta)[1][tau[75]]-widehat(beta)[1][tau[25]]))


ggsave(file="plots/S3-lineardispersion-difference-histogram.tiff", hist.coef.split, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/S3-lineardispersion-difference-histogram.jpeg", hist.coef.split, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)



## Boxplot of all diffs #####

box.coef.n <- ggplot(data=est.n.diff, aes(x=diff, y=n_obs, color=n_obs, fill=n_obs)) +
  geom_boxplot(alpha=0.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=truediff),
             size=0.5, linetype="dashed") +
  scale_color_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1) +
  scale_fill_viridis(discrete = TRUE, option="inferno", begin=0.2, end=0.7, direction=-1)


box.coef.split <- box.coef.n %+% subset(est.n.diff, !(coef.type %in% "(Intercept)")) +
  facet_wrap(~coef.type, nrow=2, ncol=2, scales="free_x") + 
  labs(x=bquote(widehat(beta)[tau[75]]-widehat(beta)[tau[25]]), y ="")

ggsave(file="plots/S3-lineardispersion-difference-boxplot.tiff", hist.coef.split, 
       device = tiff, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)

ggsave(file="plots/S3-lineardispersion-difference-boxplot.jpeg", hist.coef.split, 
       device = jpeg, width = 114.3, height = 88.9, units = "mm", limitsize = FALSE, dpi = 600)
