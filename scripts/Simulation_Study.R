##############HEADER##############################
## author:    Daniel Gotthardt
## contact:   daniel.gotthardt@studium.uni-hamburg.de / daniel.gotthardt@gmx.de
## file name: Simulation_Study.R
## Context:   Master Thesis Sociology at the Universität Hamburg
## Input:     -
## Output:    Various simulated data and additional data.frames with extracted 
##            values from the models and simulation summaries for the different szenarios
## Summary:   This R-File simulates different szenarios to evaluate models
##            for the analysis of dispersion.

## This simulation study follows broadly the advice by Morris et al 2019

################BODY##############################


# Clean up workspace -------------------------------------
rm(list=ls(all.names = TRUE))
gc()

# Setup packages once -------------------------------------
#uncomment to install the packages if necessary
#install.packages("statmod")
#install.packages("dglm")
#install.packages("quantreg")
#install.packages("rsimsum")
#install.packages("MASS)

# Packages for model fitting
library(statmod)
library(dglm)
library(quantreg)
# Packages for random data generation and simulation summaries
library(rsimsum)
library(MASS)


# Convenience Functions ------------------------------------------------
#tracking interactively if randomseed changes

if (interactive()) {
  invisible(addTaskCallback(local({
    last <- .GlobalEnv$.Random.seed
    
    function(...) {
      curr <- .GlobalEnv$.Random.seed
      if (!identical(curr, last)) {
        msg <- "NOTE: .Random.seed changed"
        if (requireNamespace("crayon", quietly=TRUE)) msg <- crayon::blurred(msg)
        message(msg)
        last <<- curr
      }
      TRUE
    }
  }), name = "RNG tracker"))
}

#define functions

## Computes a confidence interval for normal distribution; Wald-CI

CI.fun <- function(coef, se, alpha){
  lo <- coef+se*qnorm(alpha/2,0,1)   # lower bound
  hi <- coef+se*qnorm(1-(alpha/2),0,1)  # upper bound
  return(data.frame(ci_lo=lo,ci_hi=hi)) # return confidence interval as a vector
}

## Check if a (list) object is a (converged) model
is.model.fun <- function(model,legal){
  ## Check if object is a legal model object
  any(ifelse(legal %in% class(model),TRUE,FALSE))
}

#set legal values for is.model
legalmod <- c("glm","lm","dglm","rq","rqs") #change this for different methods

# Set seed once for replicability-------------------------------

set.seed(42) # change this to redo the simulation with different random numbers

# Szenario 1 Doppel-Null --------------------------------------------------------------
# Null is true for scale and mean model
## 1.1 N-obs = 30 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 30
n_obs.30 <- 30

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0
lambda_0 <- 0.5
lambda_1 <- 0
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.30 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.30 <- .Random.seed

# Draw random data
simdata.30 <- lapply(simdata.30,function(r) {
  x1 <- runif(n_obs.30, 1, 10)
  normal <- rnorm(n_obs.30,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0+x1*beta_1 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.30 <- lapply(simdata.30,function(data) data$random)
random.30 <- c(list(random0.30),random.30)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.30 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.30 <- lapply(simdata.30, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),is.model.fun,legalmod)


# Fit variance function regression with ml

ml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.30 <- !sapply(ml.fit.30,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.30 <- !sapply(reml.fit.30,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.30 <- matrix(data=NA,nrow=n_sim)
qr.fit.30 <- lapply(tau,function(t) lapply(simdata.30, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.30 <- cbind(two.step.lm.fit.na.30,two.step.gamma.fit.na.30,ml.fit.na.30,reml.fit.na.30,qr.t25.fit.na.30,qr.t50.fit.na.30,qr.t75.fit.na.30)
colSums(fit.na.30)

# 5 gamma na, 9 ml and 9 reml

# Delete simdata to save RAM

save(simdata.30, file = "dataS1_30_simdata.rds")

rm(list=c("simdata.30"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.30[,3:4] <- NA #not estimated
qr.coef.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = (z_75-z_25)*sigma=1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.30[,3:4] <- (qr.coef.q.30[,3:4]-qr.coef.q.30[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.30 <- as.data.frame(rbind(two.step.coef.30,ml.coef.30,reml.coef.30, qr.coef.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.30 <- cbind(coef.30,repetition,model)

# Extract model SE

two.step.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.30[,3:4] <- NA #not estimated
qr.se.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.30 <- as.data.frame(rbind(two.step.se.30,ml.se.30,reml.se.30,qr.se.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.30 <- cbind(se.30,repetition,model)


# Extract p-values

two.step.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.30[,3:4] <- NA #not estimated
qr.p.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.30[i,4]<- anova(qr.fit.30[["tau = 0.25"]][[i]],qr.fit.30[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.30 <- as.data.frame(rbind(two.step.p.30,ml.p.30,reml.p.30,qr.p.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.30 <- cbind(p.30,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.30 <- merge(coef.30, se.30, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.30 <- merge(est.30, p.30, by=c("repetition", "model"), sort=FALSE)
colnames(est.30)[colnames(est.30) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.30 <- reshape(est.30, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = estimands
)
# Include true values
est.long.30$true <- sapply(est.long.30$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)


# Save error fits together with repetition number

error.30 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.30 | two.step.gamma.fit.na.30),
                                  fit = two.step.fit.30[two.step.lm.fit.na.30 | two.step.gamma.fit.na.30]),
                 ML = cbind(repetition = which(ml.fit.na.30),
                            fit = ml.fit.30[ml.fit.na.30]),
                 REML = cbind(repetition = which(reml.fit.na.30),
                              fit = reml.fit.30[reml.fit.na.30]),
                 QR = cbind(repetition = which(qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30),
                            fit = qr.fit.30[qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30]))
saveRDS(error.30, file = "dataS1_30_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.30", "ml.fit.30", "reml.fit.30", "qr.fit.30"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.30$bias <- est.long.30$coef - est.long.30$true

# Coverage
# Construct CI first
conf.long.30 <- CI.fun(coef = est.long.30[,"coef"], se = est.long.30[,"se"], alpha=0.05)
est.long.30 <- cbind(est.long.30, conf.long.30)
# Check if true value is in the CI
est.long.30$cover <- (est.long.30$true >= est.long.30$ci_lo) & (est.long.30$true <= est.long.30$ci_hi)

# Rejection
reject <- est.long.30$p <= 0.05
est.long.30 <- cbind(est.long.30,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.30.mean <- aggregate(. ~ model+coef.type, data = est.long.30, na.action = na.pass,
                         FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.30, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.30.agg <- merge(est.30.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.30 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.30, na.action = na.pass, var, na.rm=T)
est.30.agg <- merge(est.30.agg, sampvar.30, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.30.agg$empse <- sqrt(est.30.agg$sampvar)


ms.30 <- multisimsum(
  data = est.long.30,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.30 <- summary(ms.30)$summ


# We need to calculate rejection MCSE on our own

est.30.agg$rejection_mcse <- sqrt(est.30.agg$reject* (1-est.30.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.30, file="S1_DoppelNull_ms_summ_30.rds")
saveRDS(est.30.agg, file="S1_DoppelNull_est_agg_30.rds")
saveRDS(est.long.30, file="S1_DoppelNull_est_long_30.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])


## 1.2 N-obs = 50 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 50
n_obs.50 <- 50

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0
lambda_0 <- 0.5
lambda_1 <- 0
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.50 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.50 <- .Random.seed

# Draw random data
simdata.50 <- lapply(simdata.50,function(r) {
  x1 <- runif(n_obs.50, 1, 10)
  normal <- rnorm(n_obs.50,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0+x1*beta_1 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.50 <- lapply(simdata.50,function(data) data$random)
random.50 <- c(list(random0.50),random.50)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.50 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.50 <- lapply(simdata.50, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.50 <- !sapply(ml.fit.50,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.50 <- !sapply(reml.fit.50,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.50 <- matrix(data=NA,nrow=n_sim)
qr.fit.50 <- lapply(tau,function(t) lapply(simdata.50, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.50 <- cbind(two.step.lm.fit.na.50,two.step.gamma.fit.na.50,ml.fit.na.50,reml.fit.na.50,qr.t25.fit.na.50,qr.t50.fit.na.50,qr.t75.fit.na.50)
colSums(fit.na.50)

# Delete simdata to save RAM

save(simdata.50, file = "dataS1_50_simdata.rds")
rm(list=c("simdata.50"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.50[,3:4] <- NA #not estimated
qr.coef.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.50[,3:4] <- (qr.coef.q.50[,3:4]-qr.coef.q.50[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.50 <- as.data.frame(rbind(two.step.coef.50,ml.coef.50,reml.coef.50, qr.coef.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.50 <- cbind(coef.50,repetition,model)

# Extract model SE

two.step.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.50[,3:4] <- NA #not estimated
qr.se.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.50 <- as.data.frame(rbind(two.step.se.50,ml.se.50,reml.se.50,qr.se.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.50 <- cbind(se.50,repetition,model)


# Extract p-values

two.step.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.50[,3:4] <- NA #not estimated
qr.p.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.50[i,4]<- anova(qr.fit.50[["tau = 0.25"]][[i]],qr.fit.50[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.50 <- as.data.frame(rbind(two.step.p.50,ml.p.50,reml.p.50,qr.p.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.50 <- cbind(p.50,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.50 <- merge(coef.50, se.50, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.50 <- merge(est.50, p.50, by=c("repetition", "model"), sort=FALSE)
colnames(est.50)[colnames(est.50) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.50 <- reshape(est.50, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = estimands
)
# Include true values
est.long.50$true <- sapply(est.long.50$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.50 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.50 | two.step.gamma.fit.na.50),
                                  fit = two.step.fit.50[two.step.lm.fit.na.50 | two.step.gamma.fit.na.50]),
                 ML = cbind(repetition = which(ml.fit.na.50),
                            fit = ml.fit.50[ml.fit.na.50]),
                 REML = cbind(repetition = which(reml.fit.na.50),
                              fit = reml.fit.50[reml.fit.na.50]),
                 QR = cbind(repetition = which(qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50),
                            fit = qr.fit.50[qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50]))
saveRDS(error.50, file = "dataS1_50_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.50", "ml.fit.50", "reml.fit.50", "qr.fit.50"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.50$bias <- est.long.50$coef - est.long.50$true

# Coverage
# Construct CI first
conf.long.50 <- CI.fun(coef = est.long.50[,"coef"], se = est.long.50[,"se"], alpha=0.05)
est.long.50 <- cbind(est.long.50, conf.long.50)
# Check if true value is in the CI
est.long.50$cover <- (est.long.50$true >= est.long.50$ci_lo) & (est.long.50$true <= est.long.50$ci_hi)

# Rejection
reject <- est.long.50$p <= 0.05
est.long.50 <- cbind(est.long.50,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.50.mean <- aggregate(. ~ model+coef.type, data = est.long.50, na.action = na.pass,
                         FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.50, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.50.agg <- merge(est.50.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.50 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.50, na.action = na.pass, var, na.rm=T)
est.50.agg <- merge(est.50.agg, sampvar.50, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.50.agg$empse <- sqrt(est.50.agg$sampvar)

ms.50 <- multisimsum(
  data = est.long.50,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.50 <- summary(ms.50)$summ


# We need to calculate rejection MCSE on our own

est.50.agg$rejection_mcse <- sqrt(est.50.agg$reject* (1-est.50.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.50, file="S1_DoppelNull_ms_summ_50.rds")
saveRDS(est.50.agg, file="S1_DoppelNull_est_agg_50.rds")
saveRDS(est.long.50, file="S1_DoppelNull_est_long_50.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

## 1.3 N-obs = 1000 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 1000
n_obs.1000 <- 1000

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0
lambda_0 <- 0.5
lambda_1 <- 0
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.1000 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.1000 <- .Random.seed

# Draw random data
simdata.1000 <- lapply(simdata.1000,function(r) {
  x1 <- runif(n_obs.1000, 1, 10)
  normal <- rnorm(n_obs.1000,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0+x1*beta_1 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.1000 <- lapply(simdata.1000,function(data) data$random)
random.1000 <- c(list(random0.1000),random.1000)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.1000 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.1000 <- lapply(simdata.1000, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.1000 <- !sapply(ml.fit.1000,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.1000 <- !sapply(reml.fit.1000,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.1000 <- matrix(data=NA,nrow=n_sim)
qr.fit.1000 <- lapply(tau,function(t) lapply(simdata.1000, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.1000 <- cbind(two.step.lm.fit.na.1000,two.step.gamma.fit.na.1000,ml.fit.na.1000,reml.fit.na.1000,qr.t25.fit.na.1000,qr.t50.fit.na.1000,qr.t75.fit.na.1000)
colSums(fit.na.1000)

# Delete simdata to save RAM

save(simdata.1000, file = "dataS1_1000_simdata.rds")
rm(list=c("simdata.1000"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.1000[,3:4] <- NA #not estimated
qr.coef.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.1000[,3:4] <- (qr.coef.q.1000[,3:4]-qr.coef.q.1000[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.1000 <- as.data.frame(rbind(two.step.coef.1000,ml.coef.1000,reml.coef.1000, qr.coef.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.1000 <- cbind(coef.1000,repetition,model)

# Extract model SE

two.step.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.1000[,3:4] <- NA #not estimated
qr.se.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.1000 <- as.data.frame(rbind(two.step.se.1000,ml.se.1000,reml.se.1000,qr.se.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.1000 <- cbind(se.1000,repetition,model)


# Extract p-values

two.step.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.1000[,3:4] <- NA #not estimated
qr.p.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.1000[i,4]<- anova(qr.fit.1000[["tau = 0.25"]][[i]],qr.fit.1000[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.1000 <- as.data.frame(rbind(two.step.p.1000,ml.p.1000,reml.p.1000,qr.p.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.1000 <- cbind(p.1000,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.1000 <- merge(coef.1000, se.1000, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.1000 <- merge(est.1000, p.1000, by=c("repetition", "model"), sort=FALSE)
colnames(est.1000)[colnames(est.1000) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.1000 <- reshape(est.1000, idvar=c("repetition","model"), direction ="long",
                         varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                         v.names = c("coef","se","p"),
                         timevar = "coef.type",
                         times = estimands
)
# Include true values
est.long.1000$true <- sapply(est.long.1000$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.1000 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000),
                                    fit = two.step.fit.1000[two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000]),
                   ML = cbind(repetition = which(ml.fit.na.1000),
                              fit = ml.fit.1000[ml.fit.na.1000]),
                   REML = cbind(repetition = which(reml.fit.na.1000),
                                fit = reml.fit.1000[reml.fit.na.1000]),
                   QR = cbind(repetition = which(qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000),
                              fit = qr.fit.1000[qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000]))
saveRDS(error.1000, file = "dataS1_1000_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.1000", "ml.fit.1000", "reml.fit.1000", "qr.fit.1000"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.1000$bias <- est.long.1000$coef - est.long.1000$true

# Coverage
# Construct CI first
conf.long.1000 <- CI.fun(coef = est.long.1000[,"coef"], se = est.long.1000[,"se"], alpha=0.05)
est.long.1000 <- cbind(est.long.1000, conf.long.1000)
# Check if true value is in the CI
est.long.1000$cover <- (est.long.1000$true >= est.long.1000$ci_lo) & (est.long.1000$true <= est.long.1000$ci_hi)

# Rejection
reject <- est.long.1000$p <= 0.05
est.long.1000 <- cbind(est.long.1000,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.1000.mean <- aggregate(. ~ model+coef.type, data = est.long.1000, na.action = na.pass,
                           FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.1000, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.1000.agg <- merge(est.1000.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.1000 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.1000, na.action = na.pass, var, na.rm=T)
est.1000.agg <- merge(est.1000.agg, sampvar.1000, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.1000.agg$empse <- sqrt(est.1000.agg$sampvar)

ms.1000 <- multisimsum(
  data = est.long.1000,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.1000 <- summary(ms.1000)$summ


# We need to calculate rejection MCSE on our own

est.1000.agg$rejection_mcse <- sqrt(est.1000.agg$reject* (1-est.1000.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.1000, file="S1_DoppelNull_ms_summ_1000.rds")
saveRDS(est.1000.agg, file="S1_DoppelNull_est_agg_1000.rds")
saveRDS(est.long.1000, file="S1_DoppelNull_est_long_1000.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

# Szenario 2 Einfach-Null --------------------------------------------------------------
# Null is true for scale model while x1 is a predictor for y1 in the mean model
## 2.1 N-obs = 30 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 30
n_obs.30 <- 30

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.5
lambda_0 <- 0.5
lambda_1 <- 0
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.30 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.30 <- .Random.seed

# Draw random data
simdata.30 <- lapply(simdata.30,function(r) {
  x1 <- runif(n_obs.30, 1, 10)
  normal <- rnorm(n_obs.30,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0+x1*beta_1 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.30 <- lapply(simdata.30,function(data) data$random)
random.30 <- c(list(random0.30),random.30)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.30 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.30 <- lapply(simdata.30, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.30 <- !sapply(ml.fit.30,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.30 <- !sapply(reml.fit.30,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.30 <- matrix(data=NA,nrow=n_sim)
qr.fit.30 <- lapply(tau,function(t) lapply(simdata.30, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.30 <- cbind(two.step.lm.fit.na.30,two.step.gamma.fit.na.30,ml.fit.na.30,reml.fit.na.30,qr.t25.fit.na.30,qr.t50.fit.na.30,qr.t75.fit.na.30)
colSums(fit.na.30)

# 8 na ML, 7 na REML, 6 na Gamma

# Delete simdata to save RAM

save(simdata.30, file = "dataS2_30_simdata.rds")
rm(list=c("simdata.30"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.30[,3:4] <- NA #not estimated
qr.coef.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.30[,3:4] <- (qr.coef.q.30[,3:4]-qr.coef.q.30[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.30 <- as.data.frame(rbind(two.step.coef.30,ml.coef.30,reml.coef.30, qr.coef.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.30 <- cbind(coef.30,repetition,model)

# Extract model SE

two.step.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.30[,3:4] <- NA #not estimated
qr.se.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.30 <- as.data.frame(rbind(two.step.se.30,ml.se.30,reml.se.30,qr.se.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.30 <- cbind(se.30,repetition,model)


# Extract p-values

two.step.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.30[,3:4] <- NA #not estimated
qr.p.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.30[i,4]<- anova(qr.fit.30[["tau = 0.25"]][[i]],qr.fit.30[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.30 <- as.data.frame(rbind(two.step.p.30,ml.p.30,reml.p.30,qr.p.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.30 <- cbind(p.30,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.30 <- merge(coef.30, se.30, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.30 <- merge(est.30, p.30, by=c("repetition", "model"), sort=FALSE)
colnames(est.30)[colnames(est.30) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.30 <- reshape(est.30, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = estimands
)
# Include true values
est.long.30$true <- sapply(est.long.30$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.30 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.30 | two.step.gamma.fit.na.30),
                                  fit = two.step.fit.30[two.step.lm.fit.na.30 | two.step.gamma.fit.na.30]),
                 ML = cbind(repetition = which(ml.fit.na.30),
                            fit = ml.fit.30[ml.fit.na.30]),
                 REML = cbind(repetition = which(reml.fit.na.30),
                              fit = reml.fit.30[reml.fit.na.30]),
                 QR = cbind(repetition = which(qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30),
                            fit = qr.fit.30[qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30]))
saveRDS(error.30, file = "dataS2_30_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.30", "ml.fit.30", "reml.fit.30", "qr.fit.30"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.30$bias <- est.long.30$coef - est.long.30$true

# Coverage
# Construct CI first
conf.long.30 <- CI.fun(coef = est.long.30[,"coef"], se = est.long.30[,"se"], alpha=0.05)
est.long.30 <- cbind(est.long.30, conf.long.30)
# Check if true value is in the CI
est.long.30$cover <- (est.long.30$true >= est.long.30$ci_lo) & (est.long.30$true <= est.long.30$ci_hi)

# Rejection
reject <- est.long.30$p <= 0.05
est.long.30 <- cbind(est.long.30,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.30.mean <- aggregate(. ~ model+coef.type, data = est.long.30, na.action = na.pass,
                         FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.30, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.30.agg <- merge(est.30.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.30 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.30, na.action = na.pass, var, na.rm=T)
est.30.agg <- merge(est.30.agg, sampvar.30, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.30.agg$empse <- sqrt(est.30.agg$sampvar)

ms.30 <- multisimsum(
  data = est.long.30,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.30 <- summary(ms.30)$summ


# We need to calculate rejection MCSE on our own

est.30.agg$rejection_mcse <- sqrt(est.30.agg$reject* (1-est.30.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.30, file="S2_Null_ms_summ_30.rds")
saveRDS(est.30.agg, file="S2_Null_est_agg_30.rds")
saveRDS(est.long.30, file="S2_Null_est_long_30.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

## 2.2 N-obs = 50 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 50
n_obs.50 <- 50

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.5
lambda_0 <- 0.5
lambda_1 <- 0
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.50 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.50 <- .Random.seed

# Draw random data
simdata.50 <- lapply(simdata.50,function(r) {
  x1 <- runif(n_obs.50, 1, 10)
  normal <- rnorm(n_obs.50,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0+x1*beta_1 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.50 <- lapply(simdata.50,function(data) data$random)
random.50 <- c(list(random0.50),random.50)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.50 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.50 <- lapply(simdata.50, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.50 <- !sapply(ml.fit.50,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.50 <- !sapply(reml.fit.50,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.50 <- matrix(data=NA,nrow=n_sim)
qr.fit.50 <- lapply(tau,function(t) lapply(simdata.50, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.50 <- cbind(two.step.lm.fit.na.50,two.step.gamma.fit.na.50,ml.fit.na.50,reml.fit.na.50,qr.t25.fit.na.50,qr.t50.fit.na.50,qr.t75.fit.na.50)
colSums(fit.na.50)


# Delete simdata to save RAM
save(simdata.50, file = "dataS2_50_simdata.rds")
rm(list=c("simdata.50"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.50[,3:4] <- NA #not estimated
qr.coef.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.50[,3:4] <- (qr.coef.q.50[,3:4]-qr.coef.q.50[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.50 <- as.data.frame(rbind(two.step.coef.50,ml.coef.50,reml.coef.50, qr.coef.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.50 <- cbind(coef.50,repetition,model)

# Extract model SE

two.step.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.50[,3:4] <- NA #not estimated
qr.se.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.50 <- as.data.frame(rbind(two.step.se.50,ml.se.50,reml.se.50,qr.se.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.50 <- cbind(se.50,repetition,model)


# Extract p-values

two.step.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.50[,3:4] <- NA #not estimated
qr.p.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.50[i,4]<- anova(qr.fit.50[["tau = 0.25"]][[i]],qr.fit.50[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.50 <- as.data.frame(rbind(two.step.p.50,ml.p.50,reml.p.50,qr.p.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.50 <- cbind(p.50,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.50 <- merge(coef.50, se.50, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.50 <- merge(est.50, p.50, by=c("repetition", "model"), sort=FALSE)
colnames(est.50)[colnames(est.50) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.50 <- reshape(est.50, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = estimands
)
# Include true values
est.long.50$true <- sapply(est.long.50$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.50 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.50 | two.step.gamma.fit.na.50),
                                  fit = two.step.fit.50[two.step.lm.fit.na.50 | two.step.gamma.fit.na.50]),
                 ML = cbind(repetition = which(ml.fit.na.50),
                            fit = ml.fit.50[ml.fit.na.50]),
                 REML = cbind(repetition = which(reml.fit.na.50),
                              fit = reml.fit.50[reml.fit.na.50]),
                 QR = cbind(repetition = which(qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50),
                            fit = qr.fit.50[qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50]))
saveRDS(error.50, file = "dataS2_50_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.50", "ml.fit.50", "reml.fit.50", "qr.fit.50"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.50$bias <- est.long.50$coef - est.long.50$true

# Coverage
# Construct CI first
conf.long.50 <- CI.fun(coef = est.long.50[,"coef"], se = est.long.50[,"se"], alpha=0.05)
est.long.50 <- cbind(est.long.50, conf.long.50)
# Check if true value is in the CI
est.long.50$cover <- (est.long.50$true >= est.long.50$ci_lo) & (est.long.50$true <= est.long.50$ci_hi)

# Rejection
reject <- est.long.50$p <= 0.05
est.long.50 <- cbind(est.long.50,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.50.mean <- aggregate(. ~ model+coef.type, data = est.long.50, na.action = na.pass,
                         FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.50, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.50.agg <- merge(est.50.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.50 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.50, na.action = na.pass, var, na.rm=T)
est.50.agg <- merge(est.50.agg, sampvar.50, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.50.agg$empse <- sqrt(est.50.agg$sampvar)

ms.50 <- multisimsum(
  data = est.long.50,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.50 <- summary(ms.50)$summ


# We need to calculate rejection MCSE on our own

est.50.agg$rejection_mcse <- sqrt(est.50.agg$reject* (1-est.50.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.50, file="S2_Null_ms_summ_50.rds")
saveRDS(est.50.agg, file="S2_Null_est_agg_50.rds")
saveRDS(est.long.50, file="S2_Null_est_long_50.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

## 2.3 N-obs = 1000 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 1000
n_obs.1000 <- 1000

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.5
lambda_0 <- 0.5
lambda_1 <- 0
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.1000 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.1000 <- .Random.seed

# Draw random data
simdata.1000 <- lapply(simdata.1000,function(r) {
  x1 <- runif(n_obs.1000, 1, 10)
  normal <- rnorm(n_obs.1000,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0+x1*beta_1 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.1000 <- lapply(simdata.1000,function(data) data$random)
random.1000 <- c(list(random0.1000),random.1000)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.1000 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.1000 <- lapply(simdata.1000, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.1000 <- !sapply(ml.fit.1000,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.1000 <- !sapply(reml.fit.1000,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.1000 <- matrix(data=NA,nrow=n_sim)
qr.fit.1000 <- lapply(tau,function(t) lapply(simdata.1000, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.1000 <- cbind(two.step.lm.fit.na.1000,two.step.gamma.fit.na.1000,ml.fit.na.1000,reml.fit.na.1000,qr.t25.fit.na.1000,qr.t50.fit.na.1000,qr.t75.fit.na.1000)
colSums(fit.na.1000)

# Delete simdata to save RAM
save(simdata.1000, file = "dataS2_1000_simdata.rds")
rm(list=c("simdata.1000"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.1000[,3:4] <- NA #not estimated
qr.coef.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.1000[,3:4] <- (qr.coef.q.1000[,3:4]-qr.coef.q.1000[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.1000 <- as.data.frame(rbind(two.step.coef.1000,ml.coef.1000,reml.coef.1000, qr.coef.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.1000 <- cbind(coef.1000,repetition,model)

# Extract model SE

two.step.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.1000[,3:4] <- NA #not estimated
qr.se.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.1000 <- as.data.frame(rbind(two.step.se.1000,ml.se.1000,reml.se.1000,qr.se.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.1000 <- cbind(se.1000,repetition,model)


# Extract p-values

two.step.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.1000[,3:4] <- NA #not estimated
qr.p.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.1000[i,4]<- anova(qr.fit.1000[["tau = 0.25"]][[i]],qr.fit.1000[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.1000 <- as.data.frame(rbind(two.step.p.1000,ml.p.1000,reml.p.1000,qr.p.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.1000 <- cbind(p.1000,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.1000 <- merge(coef.1000, se.1000, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.1000 <- merge(est.1000, p.1000, by=c("repetition", "model"), sort=FALSE)
colnames(est.1000)[colnames(est.1000) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.1000 <- reshape(est.1000, idvar=c("repetition","model"), direction ="long",
                         varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                         v.names = c("coef","se","p"),
                         timevar = "coef.type",
                         times = estimands
)
# Include true values
est.long.1000$true <- sapply(est.long.1000$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.1000 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000),
                                    fit = two.step.fit.1000[two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000]),
                   ML = cbind(repetition = which(ml.fit.na.1000),
                              fit = ml.fit.1000[ml.fit.na.1000]),
                   REML = cbind(repetition = which(reml.fit.na.1000),
                                fit = reml.fit.1000[reml.fit.na.1000]),
                   QR = cbind(repetition = which(qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000),
                              fit = qr.fit.1000[qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000]))
saveRDS(error.1000, file = "dataS2_1000_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.1000", "ml.fit.1000", "reml.fit.1000", "qr.fit.1000"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.1000$bias <- est.long.1000$coef - est.long.1000$true

# Coverage
# Construct CI first
conf.long.1000 <- CI.fun(coef = est.long.1000[,"coef"], se = est.long.1000[,"se"], alpha=0.05)
est.long.1000 <- cbind(est.long.1000, conf.long.1000)
# Check if true value is in the CI
est.long.1000$cover <- (est.long.1000$true >= est.long.1000$ci_lo) & (est.long.1000$true <= est.long.1000$ci_hi)

# Rejection
reject <- est.long.1000$p <= 0.05
est.long.1000 <- cbind(est.long.1000,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.1000.mean <- aggregate(. ~ model+coef.type, data = est.long.1000, na.action = na.pass,
                           FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.1000, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.1000.agg <- merge(est.1000.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.1000 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.1000, na.action = na.pass, var, na.rm=T)
est.1000.agg <- merge(est.1000.agg, sampvar.1000, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.1000.agg$empse <- sqrt(est.1000.agg$sampvar)

ms.1000 <- multisimsum(
  data = est.long.1000,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.1000 <- summary(ms.1000)$summ


# We need to calculate rejection MCSE on our own

est.1000.agg$rejection_mcse <- sqrt(est.1000.agg$reject* (1-est.1000.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.1000, file="S2_Null_ms_summ_1000.rds")
saveRDS(est.1000.agg, file="S2_Null_est_agg_1000.rds")
saveRDS(est.long.1000, file="S2_Null_est_long_1000.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

# SZENARIO 3 Einfach-Null aber nicht-linear im MW  --------------------------------------------------------------
# Null is true for scale model while x1 and x1^2 are predictors for y1 in the mean model
## 3.1 N-obs = 30 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 30
n_obs.30 <- 30

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.5
beta_2 <- 1
lambda_0 <- 0.5
lambda_1 <- 0
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.30 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.30 <- .Random.seed

# Draw random data
simdata.30 <- lapply(simdata.30,function(r) {
  x1 <- runif(n_obs.30, 1, 10)
  normal <- rnorm(n_obs.30,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0+x1*beta_1+x1^2*beta_2 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.30 <- lapply(simdata.30,function(data) data$random)
random.30 <- c(list(random0.30),random.30)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.30 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.30 <- lapply(simdata.30, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.30 <- !sapply(ml.fit.30,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.30 <- !sapply(reml.fit.30,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.30 <- matrix(data=NA,nrow=n_sim)
qr.fit.30 <- lapply(tau,function(t) lapply(simdata.30, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.30 <- cbind(two.step.lm.fit.na.30,two.step.gamma.fit.na.30,ml.fit.na.30,reml.fit.na.30,qr.t25.fit.na.30,qr.t50.fit.na.30,qr.t75.fit.na.30)
colSums(fit.na.30)

# 95 ML, 89 REML, 106 Gamma

# Delete simdata to save RAM
save(simdata.30, file = "dataS3_30_simdata.rds")
rm(list=c("simdata.30"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.30[,3:4] <- NA #not estimated
qr.coef.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.30[,3:4] <- (qr.coef.q.30[,3:4]-qr.coef.q.30[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.30 <- as.data.frame(rbind(two.step.coef.30,ml.coef.30,reml.coef.30, qr.coef.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.30 <- cbind(coef.30,repetition,model)

# Extract model SE

two.step.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.30[,3:4] <- NA #not estimated
qr.se.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.30 <- as.data.frame(rbind(two.step.se.30,ml.se.30,reml.se.30,qr.se.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.30 <- cbind(se.30,repetition,model)


# Extract p-values

two.step.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.30[,3:4] <- NA #not estimated
qr.p.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.30[i,4]<- anova(qr.fit.30[["tau = 0.25"]][[i]],qr.fit.30[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.30 <- as.data.frame(rbind(two.step.p.30,ml.p.30,reml.p.30,qr.p.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.30 <- cbind(p.30,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.30 <- merge(coef.30, se.30, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.30 <- merge(est.30, p.30, by=c("repetition", "model"), sort=FALSE)
colnames(est.30)[colnames(est.30) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.30 <- reshape(est.30, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = estimands
)
# Include true values
est.long.30$true <- sapply(est.long.30$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.30 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.30 | two.step.gamma.fit.na.30),
                                  fit = two.step.fit.30[two.step.lm.fit.na.30 | two.step.gamma.fit.na.30]),
                 ML = cbind(repetition = which(ml.fit.na.30),
                            fit = ml.fit.30[ml.fit.na.30]),
                 REML = cbind(repetition = which(reml.fit.na.30),
                              fit = reml.fit.30[reml.fit.na.30]),
                 QR = cbind(repetition = which(qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30),
                            fit = qr.fit.30[qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30]))
saveRDS(error.30, file = "dataS3_30_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.30", "ml.fit.30", "reml.fit.30", "qr.fit.30"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.30$bias <- est.long.30$coef - est.long.30$true

# Coverage
# Construct CI first
conf.long.30 <- CI.fun(coef = est.long.30[,"coef"], se = est.long.30[,"se"], alpha=0.05)
est.long.30 <- cbind(est.long.30, conf.long.30)
# Check if true value is in the CI
est.long.30$cover <- (est.long.30$true >= est.long.30$ci_lo) & (est.long.30$true <= est.long.30$ci_hi)

# Rejection
reject <- est.long.30$p <= 0.05
est.long.30 <- cbind(est.long.30,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.30.mean <- aggregate(. ~ model+coef.type, data = est.long.30, na.action = na.pass,
                         FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.30, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.30.agg <- merge(est.30.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.30 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.30, na.action = na.pass, var, na.rm=T)
est.30.agg <- merge(est.30.agg, sampvar.30, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.30.agg$empse <- sqrt(est.30.agg$sampvar)

ms.30 <- multisimsum(
  data = est.long.30,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.30 <- summary(ms.30)$summ


# We need to calculate rejection MCSE on our own

est.30.agg$rejection_mcse <- sqrt(est.30.agg$reject* (1-est.30.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.30, file="S3_Nullx2_ms_summ_30.rds")
saveRDS(est.30.agg, file="S3_Nullx2_est_agg_30.rds")
saveRDS(est.long.30, file="S3_Nullx2_est_long_30.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

## 3.2 N-obs = 50 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 50
n_obs.50 <- 50

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.5
beta_2 <- 1
lambda_0 <- 0.5
lambda_1 <- 0
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.50 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.50 <- .Random.seed

# Draw random data
simdata.50 <- lapply(simdata.50,function(r) {
  x1 <- runif(n_obs.50, 1, 10)
  normal <- rnorm(n_obs.50,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0+x1*beta_1+x1^2*beta_2 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.50 <- lapply(simdata.50,function(data) data$random)
random.50 <- c(list(random0.50),random.50)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.50 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.50 <- lapply(simdata.50, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.50 <- !sapply(ml.fit.50,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.50 <- !sapply(reml.fit.50,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.50 <- matrix(data=NA,nrow=n_sim)
qr.fit.50 <- lapply(tau,function(t) lapply(simdata.50, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.50 <- cbind(two.step.lm.fit.na.50,two.step.gamma.fit.na.50,ml.fit.na.50,reml.fit.na.50,qr.t25.fit.na.50,qr.t50.fit.na.50,qr.t75.fit.na.50)
colSums(fit.na.50)

# 68 ML, 85 REML, 37 Gamma


# Delete simdata to save RAM
save(simdata.50, file = "dataS3_50_simdata.rds")
rm(list=c("simdata.50"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.50[,3:4] <- NA #not estimated
qr.coef.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.50[,3:4] <- (qr.coef.q.50[,3:4]-qr.coef.q.50[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.50 <- as.data.frame(rbind(two.step.coef.50,ml.coef.50,reml.coef.50, qr.coef.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.50 <- cbind(coef.50,repetition,model)

# Extract model SE

two.step.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.50[,3:4] <- NA #not estimated
qr.se.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.50 <- as.data.frame(rbind(two.step.se.50,ml.se.50,reml.se.50,qr.se.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.50 <- cbind(se.50,repetition,model)


# Extract p-values

two.step.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.50[,3:4] <- NA #not estimated
qr.p.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.50[i,4]<- anova(qr.fit.50[["tau = 0.25"]][[i]],qr.fit.50[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.50 <- as.data.frame(rbind(two.step.p.50,ml.p.50,reml.p.50,qr.p.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.50 <- cbind(p.50,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.50 <- merge(coef.50, se.50, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.50 <- merge(est.50, p.50, by=c("repetition", "model"), sort=FALSE)
colnames(est.50)[colnames(est.50) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.50 <- reshape(est.50, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = estimands
)
# Include true values
est.long.50$true <- sapply(est.long.50$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.50 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.50 | two.step.gamma.fit.na.50),
                                  fit = two.step.fit.50[two.step.lm.fit.na.50 | two.step.gamma.fit.na.50]),
                 ML = cbind(repetition = which(ml.fit.na.50),
                            fit = ml.fit.50[ml.fit.na.50]),
                 REML = cbind(repetition = which(reml.fit.na.50),
                              fit = reml.fit.50[reml.fit.na.50]),
                 QR = cbind(repetition = which(qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50),
                            fit = qr.fit.50[qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50]))
saveRDS(error.50, file = "dataS3_50_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.50", "ml.fit.50", "reml.fit.50", "qr.fit.50"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.50$bias <- est.long.50$coef - est.long.50$true

# Coverage
# Construct CI first
conf.long.50 <- CI.fun(coef = est.long.50[,"coef"], se = est.long.50[,"se"], alpha=0.05)
est.long.50 <- cbind(est.long.50, conf.long.50)
# Check if true value is in the CI
est.long.50$cover <- (est.long.50$true >= est.long.50$ci_lo) & (est.long.50$true <= est.long.50$ci_hi)

# Rejection
reject <- est.long.50$p <= 0.05
est.long.50 <- cbind(est.long.50,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.50.mean <- aggregate(. ~ model+coef.type, data = est.long.50, na.action = na.pass,
                         FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.50, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.50.agg <- merge(est.50.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.50 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.50, na.action = na.pass, var, na.rm=T)
est.50.agg <- merge(est.50.agg, sampvar.50, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.50.agg$empse <- sqrt(est.50.agg$sampvar)

ms.50 <- multisimsum(
  data = est.long.50,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.50 <- summary(ms.50)$summ


# We need to calculate rejection MCSE on our own

est.50.agg$rejection_mcse <- sqrt(est.50.agg$reject* (1-est.50.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.50, file="S3_Nullx2_ms_summ_50.rds")
saveRDS(est.50.agg, file="S3_Nullx2_est_agg_50.rds")
saveRDS(est.long.50, file="S3_Nullx2_est_long_50.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

## 3.3 N-obs = 1000 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 1000
n_obs.1000 <- 1000

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.5
beta_2 <- 1
lambda_0 <- 0.5
lambda_1 <- 0
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.1000 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.1000 <- .Random.seed

# Draw random data
simdata.1000 <- lapply(simdata.1000,function(r) {
  x1 <- runif(n_obs.1000, 1, 10)
  normal <- rnorm(n_obs.1000,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0+x1*beta_1+x1^2*beta_2 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.1000 <- lapply(simdata.1000,function(data) data$random)
random.1000 <- c(list(random0.1000),random.1000)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.1000 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.1000 <- lapply(simdata.1000, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.1000 <- !sapply(ml.fit.1000,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.1000 <- !sapply(reml.fit.1000,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.1000 <- matrix(data=NA,nrow=n_sim)
qr.fit.1000 <- lapply(tau,function(t) lapply(simdata.1000, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.1000 <- cbind(two.step.lm.fit.na.1000,two.step.gamma.fit.na.1000,ml.fit.na.1000,reml.fit.na.1000,qr.t25.fit.na.1000,qr.t50.fit.na.1000,qr.t75.fit.na.1000)
colSums(fit.na.1000)

# Delete simdata to save RAM
save(simdata.1000, file = "dataS3_1000_simdata.rds")
rm(list=c("simdata.1000"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.1000[,3:4] <- NA #not estimated
qr.coef.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.1000[,3:4] <- (qr.coef.q.1000[,3:4]-qr.coef.q.1000[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.1000 <- as.data.frame(rbind(two.step.coef.1000,ml.coef.1000,reml.coef.1000, qr.coef.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.1000 <- cbind(coef.1000,repetition,model)

# Extract model SE

two.step.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.1000[,3:4] <- NA #not estimated
qr.se.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.1000 <- as.data.frame(rbind(two.step.se.1000,ml.se.1000,reml.se.1000,qr.se.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.1000 <- cbind(se.1000,repetition,model)


# Extract p-values

two.step.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.1000[,3:4] <- NA #not estimated
qr.p.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.1000[i,4]<- anova(qr.fit.1000[["tau = 0.25"]][[i]],qr.fit.1000[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.1000 <- as.data.frame(rbind(two.step.p.1000,ml.p.1000,reml.p.1000,qr.p.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.1000 <- cbind(p.1000,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.1000 <- merge(coef.1000, se.1000, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.1000 <- merge(est.1000, p.1000, by=c("repetition", "model"), sort=FALSE)
colnames(est.1000)[colnames(est.1000) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.1000 <- reshape(est.1000, idvar=c("repetition","model"), direction ="long",
                         varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                         v.names = c("coef","se","p"),
                         timevar = "coef.type",
                         times = estimands
)
# Include true values
est.long.1000$true <- sapply(est.long.1000$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.1000 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000),
                                    fit = two.step.fit.1000[two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000]),
                   ML = cbind(repetition = which(ml.fit.na.1000),
                              fit = ml.fit.1000[ml.fit.na.1000]),
                   REML = cbind(repetition = which(reml.fit.na.1000),
                                fit = reml.fit.1000[reml.fit.na.1000]),
                   QR = cbind(repetition = which(qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000),
                              fit = qr.fit.1000[qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000]))
saveRDS(error.1000, file = "dataS3_1000_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.1000", "ml.fit.1000", "reml.fit.1000", "qr.fit.1000"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.1000$bias <- est.long.1000$coef - est.long.1000$true

# Coverage
# Construct CI first
conf.long.1000 <- CI.fun(coef = est.long.1000[,"coef"], se = est.long.1000[,"se"], alpha=0.05)
est.long.1000 <- cbind(est.long.1000, conf.long.1000)
# Check if true value is in the CI
est.long.1000$cover <- (est.long.1000$true >= est.long.1000$ci_lo) & (est.long.1000$true <= est.long.1000$ci_hi)

# Rejection
reject <- est.long.1000$p <= 0.05
est.long.1000 <- cbind(est.long.1000,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.1000.mean <- aggregate(. ~ model+coef.type, data = est.long.1000, na.action = na.pass,
                           FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.1000, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.1000.agg <- merge(est.1000.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.1000 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.1000, na.action = na.pass, var, na.rm=T)
est.1000.agg <- merge(est.1000.agg, sampvar.1000, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.1000.agg$empse <- sqrt(est.1000.agg$sampvar)

ms.1000 <- multisimsum(
  data = est.long.1000,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.1000 <- summary(ms.1000)$summ


# We need to calculate rejection MCSE on our own

est.1000.agg$rejection_mcse <- sqrt(est.1000.agg$reject* (1-est.1000.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.1000, file="S3_Nullx2_ms_summ_1000.rds")
saveRDS(est.1000.agg, file="S3_Nullx2_est_agg_1000.rds")
saveRDS(est.long.1000, file="S3_Nullx2_est_long_1000.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])


# SZENARIO 4 Double Null with vertical outliers  --------------------------------------------------------------
# Null is true for scale and mean model but random vertical outliers contaminate the data
## 4.1 N-obs = 30 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 30
n_obs.30 <- 30

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.0
lambda_0 <- 0.5
lambda_1 <- 0
k <- 0.1 # 10% outlier
o <- 3 # inflation of sd for outlier
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.30 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.30 <- .Random.seed

# Draw random data
simdata.30 <- lapply(simdata.30,function(r) {
  x1 <- runif(n_obs.30, 1, 10)
  normal <- rnorm(n_obs.30*(1-k), mean=0, sd=sigma*o*sqrt(k/(1-k))) #reduced variance, so that the variance of the mixture is still sigma
  outlier <- rnorm(n_obs.30*k, mean=0,sd=3*sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*c(normal,outlier)
  y1 <- beta_0+x1*beta_1+ u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.30 <- lapply(simdata.30,function(data) data$random)
random.30 <- c(list(random0.30),random.30)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")


### Fit regressions #########################

# Fit two-step regression

two.step.fit.30 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.30 <- lapply(simdata.30, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.30 <- !sapply(ml.fit.30,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.30 <- !sapply(reml.fit.30,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.30 <- matrix(data=NA,nrow=n_sim)
qr.fit.30 <- lapply(tau,function(t) lapply(simdata.30, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.30 <- cbind(two.step.lm.fit.na.30,two.step.gamma.fit.na.30,ml.fit.na.30,reml.fit.na.30,qr.t25.fit.na.30,qr.t50.fit.na.30,qr.t75.fit.na.30)
colSums(fit.na.30)

# Delete simdata to save RAM
save(simdata.30, file = "dataS4_30_simdata.rds")
rm(list=c("simdata.30"))

# 23 na ML, 18 na REML, 18 na GAMMA


### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.30[,3:4] <- NA #not estimated
qr.coef.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.30[,3:4] <- (qr.coef.q.30[,3:4]-qr.coef.q.30[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.30 <- as.data.frame(rbind(two.step.coef.30,ml.coef.30,reml.coef.30, qr.coef.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.30 <- cbind(coef.30,repetition,model)

# Extract model SE

two.step.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.30[,3:4] <- NA #not estimated
qr.se.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.30 <- as.data.frame(rbind(two.step.se.30,ml.se.30,reml.se.30,qr.se.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.30 <- cbind(se.30,repetition,model)


# Extract p-values

two.step.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.30[,3:4] <- NA #not estimated
qr.p.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.30[i,4]<- anova(qr.fit.30[["tau = 0.25"]][[i]],qr.fit.30[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.30 <- as.data.frame(rbind(two.step.p.30,ml.p.30,reml.p.30,qr.p.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.30 <- cbind(p.30,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.30 <- merge(coef.30, se.30, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.30 <- merge(est.30, p.30, by=c("repetition", "model"), sort=FALSE)
colnames(est.30)[colnames(est.30) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.30 <- reshape(est.30, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = estimands
)
# Include true values
est.long.30$true <- sapply(est.long.30$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.30 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.30 | two.step.gamma.fit.na.30),
                                  fit = two.step.fit.30[two.step.lm.fit.na.30 | two.step.gamma.fit.na.30]),
                 ML = cbind(repetition = which(ml.fit.na.30),
                            fit = ml.fit.30[ml.fit.na.30]),
                 REML = cbind(repetition = which(reml.fit.na.30),
                              fit = reml.fit.30[reml.fit.na.30]),
                 QR = cbind(repetition = which(qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30),
                            fit = qr.fit.30[qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30]))
saveRDS(error.30, file = "dataS4_30_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.30", "ml.fit.30", "reml.fit.30", "qr.fit.30"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.30$bias <- est.long.30$coef - est.long.30$true

# Coverage
# Construct CI first
conf.long.30 <- CI.fun(coef = est.long.30[,"coef"], se = est.long.30[,"se"], alpha=0.05)
est.long.30 <- cbind(est.long.30, conf.long.30)
# Check if true value is in the CI
est.long.30$cover <- (est.long.30$true >= est.long.30$ci_lo) & (est.long.30$true <= est.long.30$ci_hi)

# Rejection
reject <- est.long.30$p <= 0.05
est.long.30 <- cbind(est.long.30,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.30.mean <- aggregate(. ~ model+coef.type, data = est.long.30, na.action = na.pass,
                         FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.30, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.30.agg <- merge(est.30.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.30 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.30, na.action = na.pass, var, na.rm=T)
est.30.agg <- merge(est.30.agg, sampvar.30, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.30.agg$empse <- sqrt(est.30.agg$sampvar)

ms.30 <- multisimsum(
  data = est.long.30,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.30 <- summary(ms.30)$summ


# We need to calculate rejection MCSE on our own

est.30.agg$rejection_mcse <- sqrt(est.30.agg$reject* (1-est.30.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.30, file="S4_Nullout_ms_summ_30.rds")
saveRDS(est.30.agg, file="S4_Nullout_est_agg_30.rds")
saveRDS(est.long.30, file="S4_Nullout_est_long_30.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])


## 4.2 N-obs = 50 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 50
n_obs.50 <- 50

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.0
lambda_0 <- 0.5
lambda_1 <- 0
k <- 0.1 # 10% outlier
o <- 3 # inflation of sd for outlier
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.50 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.50 <- .Random.seed

# Draw random data
simdata.50 <- lapply(simdata.50,function(r) {
  x1 <- runif(n_obs.50, 1, 10)
  normal <- rnorm(n_obs.50*(1-k), mean=0, sd=sigma*o*sqrt(k/(1-k))) #reduced variance, so that the variance of the mixture is still sigma
  outlier <- rnorm(n_obs.50*k, mean=0,sd=3*sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*c(normal,outlier)
  y1 <- beta_0+x1*beta_1+ u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.50 <- lapply(simdata.50,function(data) data$random)
random.50 <- c(list(random0.50),random.50)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.50 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.50 <- lapply(simdata.50, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.50 <- !sapply(ml.fit.50,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.50 <- !sapply(reml.fit.50,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.50 <- matrix(data=NA,nrow=n_sim)
qr.fit.50 <- lapply(tau,function(t) lapply(simdata.50, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.50 <- cbind(two.step.lm.fit.na.50,two.step.gamma.fit.na.50,ml.fit.na.50,reml.fit.na.50,qr.t25.fit.na.50,qr.t50.fit.na.50,qr.t75.fit.na.50)
colSums(fit.na.50)

# 11 na ML, 11 na REML, gamma 9

# Delete simdata to save RAM
save(simdata.50, file = "dataS4_50_simdata.rds")
rm(list=c("simdata.50"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.50[,3:4] <- NA #not estimated
qr.coef.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.50[,3:4] <- (qr.coef.q.50[,3:4]-qr.coef.q.50[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.50 <- as.data.frame(rbind(two.step.coef.50,ml.coef.50,reml.coef.50, qr.coef.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.50 <- cbind(coef.50,repetition,model)

# Extract model SE

two.step.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.50[,3:4] <- NA #not estimated
qr.se.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.50 <- as.data.frame(rbind(two.step.se.50,ml.se.50,reml.se.50,qr.se.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.50 <- cbind(se.50,repetition,model)


# Extract p-values

two.step.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.50[,3:4] <- NA #not estimated
qr.p.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.50[i,4]<- anova(qr.fit.50[["tau = 0.25"]][[i]],qr.fit.50[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.50 <- as.data.frame(rbind(two.step.p.50,ml.p.50,reml.p.50,qr.p.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.50 <- cbind(p.50,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.50 <- merge(coef.50, se.50, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.50 <- merge(est.50, p.50, by=c("repetition", "model"), sort=FALSE)
colnames(est.50)[colnames(est.50) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.50 <- reshape(est.50, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = estimands
)
# Include true values
est.long.50$true <- sapply(est.long.50$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.50 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.50 | two.step.gamma.fit.na.50),
                                  fit = two.step.fit.50[two.step.lm.fit.na.50 | two.step.gamma.fit.na.50]),
                 ML = cbind(repetition = which(ml.fit.na.50),
                            fit = ml.fit.50[ml.fit.na.50]),
                 REML = cbind(repetition = which(reml.fit.na.50),
                              fit = reml.fit.50[reml.fit.na.50]),
                 QR = cbind(repetition = which(qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50),
                            fit = qr.fit.50[qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50]))
saveRDS(error.50, file = "dataS4_50_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.50", "ml.fit.50", "reml.fit.50", "qr.fit.50"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.50$bias <- est.long.50$coef - est.long.50$true

# Coverage
# Construct CI first
conf.long.50 <- CI.fun(coef = est.long.50[,"coef"], se = est.long.50[,"se"], alpha=0.05)
est.long.50 <- cbind(est.long.50, conf.long.50)
# Check if true value is in the CI
est.long.50$cover <- (est.long.50$true >= est.long.50$ci_lo) & (est.long.50$true <= est.long.50$ci_hi)

# Rejection
reject <- est.long.50$p <= 0.05
est.long.50 <- cbind(est.long.50,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.50.mean <- aggregate(. ~ model+coef.type, data = est.long.50, na.action = na.pass,
                         FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.50, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.50.agg <- merge(est.50.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.50 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.50, na.action = na.pass, var, na.rm=T)
est.50.agg <- merge(est.50.agg, sampvar.50, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.50.agg$empse <- sqrt(est.50.agg$sampvar)

ms.50 <- multisimsum(
  data = est.long.50,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.50 <- summary(ms.50)$summ


# We need to calculate rejection MCSE on our own

est.50.agg$rejection_mcse <- sqrt(est.50.agg$reject* (1-est.50.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.50, file="S4_Nullout_ms_summ_50.rds")
saveRDS(est.50.agg, file="S4_Nullout_est_agg_50.rds")
saveRDS(est.long.50, file="S4_Nullout_est_long_50.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])


## 4.3 N-obs = 1000 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 1000
n_obs.1000 <- 1000

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.0
lambda_0 <- 0.5
lambda_1 <- 0
k <- 0.1 # 10% outlier
o <- 3 # inflation of sd for outlier
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.1000 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.1000 <- .Random.seed

# Draw random data
simdata.1000 <- lapply(simdata.1000,function(r) {
  x1 <- runif(n_obs.1000, 1, 10)
  normal <- rnorm(n_obs.1000*(1-k), mean=0, sd=sigma*o*sqrt(k/(1-k))) #reduced variance, so that the variance of the mixture is still sigma
  outlier <- rnorm(n_obs.1000*k, mean=0,sd=3*sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*c(normal,outlier)
  y1 <- beta_0+x1*beta_1+ u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.1000 <- lapply(simdata.1000,function(data) data$random)
random.1000 <- c(list(random0.1000),random.1000)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.1000 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.1000 <- lapply(simdata.1000, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.1000 <- !sapply(ml.fit.1000,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.1000 <- !sapply(reml.fit.1000,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.1000 <- matrix(data=NA,nrow=n_sim)
qr.fit.1000 <- lapply(tau,function(t) lapply(simdata.1000, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.1000 <- cbind(two.step.lm.fit.na.1000,two.step.gamma.fit.na.1000,ml.fit.na.1000,reml.fit.na.1000,qr.t25.fit.na.1000,qr.t50.fit.na.1000,qr.t75.fit.na.1000)
colSums(fit.na.1000)



# Delete simdata to save RAM
save(simdata.1000, file = "dataS4_1000_simdata.rds")
rm(list=c("simdata.1000"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.1000[,3:4] <- NA #not estimated
qr.coef.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.1000[,3:4] <- (qr.coef.q.1000[,3:4]-qr.coef.q.1000[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.1000 <- as.data.frame(rbind(two.step.coef.1000,ml.coef.1000,reml.coef.1000, qr.coef.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.1000 <- cbind(coef.1000,repetition,model)

# Extract model SE

two.step.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.1000[,3:4] <- NA #not estimated
qr.se.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.1000 <- as.data.frame(rbind(two.step.se.1000,ml.se.1000,reml.se.1000,qr.se.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.1000 <- cbind(se.1000,repetition,model)


# Extract p-values

two.step.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.1000[,3:4] <- NA #not estimated
qr.p.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.1000[i,4]<- anova(qr.fit.1000[["tau = 0.25"]][[i]],qr.fit.1000[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.1000 <- as.data.frame(rbind(two.step.p.1000,ml.p.1000,reml.p.1000,qr.p.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.1000 <- cbind(p.1000,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.1000 <- merge(coef.1000, se.1000, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.1000 <- merge(est.1000, p.1000, by=c("repetition", "model"), sort=FALSE)
colnames(est.1000)[colnames(est.1000) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.1000 <- reshape(est.1000, idvar=c("repetition","model"), direction ="long",
                         varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                         v.names = c("coef","se","p"),
                         timevar = "coef.type",
                         times = estimands
)
# Include true values
est.long.1000$true <- sapply(est.long.1000$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.1000 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000),
                                    fit = two.step.fit.1000[two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000]),
                   ML = cbind(repetition = which(ml.fit.na.1000),
                              fit = ml.fit.1000[ml.fit.na.1000]),
                   REML = cbind(repetition = which(reml.fit.na.1000),
                                fit = reml.fit.1000[reml.fit.na.1000]),
                   QR = cbind(repetition = which(qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000),
                              fit = qr.fit.1000[qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000]))
saveRDS(error.1000, file = "dataS4_1000_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.1000", "ml.fit.1000", "reml.fit.1000", "qr.fit.1000"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.1000$bias <- est.long.1000$coef - est.long.1000$true

# Coverage
# Construct CI first
conf.long.1000 <- CI.fun(coef = est.long.1000[,"coef"], se = est.long.1000[,"se"], alpha=0.05)
est.long.1000 <- cbind(est.long.1000, conf.long.1000)
# Check if true value is in the CI
est.long.1000$cover <- (est.long.1000$true >= est.long.1000$ci_lo) & (est.long.1000$true <= est.long.1000$ci_hi)

# Rejection
reject <- est.long.1000$p <= 0.05
est.long.1000 <- cbind(est.long.1000,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.1000.mean <- aggregate(. ~ model+coef.type, data = est.long.1000, na.action = na.pass,
                           FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.1000, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.1000.agg <- merge(est.1000.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.1000 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.1000, na.action = na.pass, var, na.rm=T)
est.1000.agg <- merge(est.1000.agg, sampvar.1000, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.1000.agg$empse <- sqrt(est.1000.agg$sampvar)

ms.1000 <- multisimsum(
  data = est.long.1000,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.1000 <- summary(ms.1000)$summ


# We need to calculate rejection MCSE on our own

est.1000.agg$rejection_mcse <- sqrt(est.1000.agg$reject* (1-est.1000.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.1000, file="S4_Nullout_ms_summ_1000.rds")
saveRDS(est.1000.agg, file="S4_Nullout_est_agg_1000.rds")
saveRDS(est.long.1000, file="S4_Nullout_est_long_1000.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])


# SZENARIO 5 Double Null with influential leverage points  --------------------------------------------------------------
# Null is true for scale and mean model but a random mass of leverage points contaminate the data
## 5.1 N-obs = 30 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 30
n_obs.30 <- 30

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.0
lambda_0 <- 0.5
lambda_1 <- 0
k <- 0.1 # 10% outlier
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.30 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.30 <- .Random.seed

# Draw random data
simdata.30 <- lapply(simdata.30,function(r) {
  x1 <- runif(n_obs.30*(1-k), 1, 10)
  normal <- rnorm(n_obs.30*(1-k), mean = 0, sd = sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0 + x1*beta_1 + u
  outlier <- mvrnorm(n_obs.30*k, mu = c(8,beta_0+3.5+8*beta_1), Sigma = matrix(c(0.5,0,0,0.5),2,2))
  x1 <- c(x1,outlier[,1])
  y1 <- c(y1,outlier[,2])
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.30 <- lapply(simdata.30,function(data) data$random)
random.30 <- c(list(random0.30),random.30)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.30 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.30 <- lapply(simdata.30, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.30 <- !sapply(ml.fit.30,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.30 <- !sapply(reml.fit.30,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.30 <- matrix(data=NA,nrow=n_sim)
qr.fit.30 <- lapply(tau,function(t) lapply(simdata.30, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.30 <- !sapply(qr.fit.30[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.30 <- cbind(two.step.lm.fit.na.30,two.step.gamma.fit.na.30,ml.fit.na.30,reml.fit.na.30,qr.t25.fit.na.30,qr.t50.fit.na.30,qr.t75.fit.na.30)
colSums(fit.na.30)

# 3 ML, 4 REML, 5 Gamma

# Delete simdata to save RAM
save(simdata.30, file = "dataS5_30_simdata.rds")
rm(list=c("simdata.30"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.30[,3:4] <- NA #not estimated
qr.coef.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.30[,3:4] <- (qr.coef.q.30[,3:4]-qr.coef.q.30[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.30 <- as.data.frame(rbind(two.step.coef.30,ml.coef.30,reml.coef.30, qr.coef.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.30 <- cbind(coef.30,repetition,model)

# Extract model SE

two.step.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.30[,3:4] <- NA #not estimated
qr.se.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.30 <- as.data.frame(rbind(two.step.se.30,ml.se.30,reml.se.30,qr.se.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.30 <- cbind(se.30,repetition,model)


# Extract p-values

two.step.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.30[,1:2] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.30[,3:4] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.30[,1:2] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.30[,3:4] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.30[,1:2] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.30[,3:4] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.30) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.30[,3:4] <- NA #not estimated
qr.p.q.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.30) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.30[,1:2] <- t(sapply(qr.fit.30[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.30[,3:4] <- t(sapply(qr.fit.30[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.30[i,4]<- anova(qr.fit.30[["tau = 0.25"]][[i]],qr.fit.30[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.30 <- as.data.frame(rbind(two.step.p.30,ml.p.30,reml.p.30,qr.p.30))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.30 <- cbind(p.30,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.30 <- merge(coef.30, se.30, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.30 <- merge(est.30, p.30, by=c("repetition", "model"), sort=FALSE)
colnames(est.30)[colnames(est.30) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.30 <- reshape(est.30, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = estimands
)
# Include true values
est.long.30$true <- sapply(est.long.30$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.30 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.30 | two.step.gamma.fit.na.30),
                                  fit = two.step.fit.30[two.step.lm.fit.na.30 | two.step.gamma.fit.na.30]),
                 ML = cbind(repetition = which(ml.fit.na.30),
                            fit = ml.fit.30[ml.fit.na.30]),
                 REML = cbind(repetition = which(reml.fit.na.30),
                              fit = reml.fit.30[reml.fit.na.30]),
                 QR = cbind(repetition = which(qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30),
                            fit = qr.fit.30[qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30]))
saveRDS(error.30, file = "dataS5_30_errors.rds")


# Delete all fits to save RAM

rm(list=c("two.step.fit.30", "ml.fit.30", "reml.fit.30", "qr.fit.30"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.30$bias <- est.long.30$coef - est.long.30$true

# Coverage
# Construct CI first
conf.long.30 <- CI.fun(coef = est.long.30[,"coef"], se = est.long.30[,"se"], alpha=0.05)
est.long.30 <- cbind(est.long.30, conf.long.30)
# Check if true value is in the CI
est.long.30$cover <- (est.long.30$true >= est.long.30$ci_lo) & (est.long.30$true <= est.long.30$ci_hi)

# Rejection
reject <- est.long.30$p <= 0.05
est.long.30 <- cbind(est.long.30,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.30.mean <- aggregate(. ~ model+coef.type, data = est.long.30, na.action = na.pass,
                         FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.30, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.30.agg <- merge(est.30.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.30 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.30, na.action = na.pass, var, na.rm=T)
est.30.agg <- merge(est.30.agg, sampvar.30, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.30.agg$empse <- sqrt(est.30.agg$sampvar)

ms.30 <- multisimsum(
  data = est.long.30,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.30 <- summary(ms.30)$summ


# We need to calculate rejection MCSE on our own

est.30.agg$rejection_mcse <- sqrt(est.30.agg$reject* (1-est.30.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.30, file="S5_Nullinf_ms_summ_30.rds")
saveRDS(est.30.agg, file="S5_Nullinf_est_agg_30.rds")
saveRDS(est.long.30, file="S5_Nullinf_est_long_30.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

## 5.2 N-obs = 50 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 50
n_obs.50 <- 50

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.0
lambda_0 <- 0.5
lambda_1 <- 0
k <- 0.1 # 10% outlier
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.50 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.50 <- .Random.seed

# Draw random data
simdata.50 <- lapply(simdata.50,function(r) {
  x1 <- runif(n_obs.50*(1-k), 1, 10)
  normal <- rnorm(n_obs.50*(1-k), mean = 0, sd = sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0 + x1*beta_1 + u
  outlier <- mvrnorm(n_obs.50*k, mu = c(8,beta_0+3.5+8*beta_1), Sigma = matrix(c(0.5,0,0,0.5),2,2))
  x1 <- c(x1,outlier[,1])
  y1 <- c(y1,outlier[,2])
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.50 <- lapply(simdata.50,function(data) data$random)
random.50 <- c(list(random0.50),random.50)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.50 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.50 <- lapply(simdata.50, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.50 <- !sapply(ml.fit.50,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.50 <- !sapply(reml.fit.50,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.50 <- matrix(data=NA,nrow=n_sim)
qr.fit.50 <- lapply(tau,function(t) lapply(simdata.50, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.50 <- !sapply(qr.fit.50[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.50 <- cbind(two.step.lm.fit.na.50,two.step.gamma.fit.na.50,ml.fit.na.50,reml.fit.na.50,qr.t25.fit.na.50,qr.t50.fit.na.50,qr.t75.fit.na.50)
colSums(fit.na.50)

# 2 na ML, 1 na REML, 1 na Two-Step

# Delete simdata to save RAM
save(simdata.50, file = "dataS5_50_simdata.rds")
rm(list=c("simdata.50"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.50[,3:4] <- NA #not estimated
qr.coef.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.50[,3:4] <- (qr.coef.q.50[,3:4]-qr.coef.q.50[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.50 <- as.data.frame(rbind(two.step.coef.50,ml.coef.50,reml.coef.50, qr.coef.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.50 <- cbind(coef.50,repetition,model)

# Extract model SE

two.step.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.50[,3:4] <- NA #not estimated
qr.se.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.50 <- as.data.frame(rbind(two.step.se.50,ml.se.50,reml.se.50,qr.se.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.50 <- cbind(se.50,repetition,model)


# Extract p-values

two.step.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.50[,1:2] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.50[,3:4] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.50[,1:2] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.50[,3:4] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.50[,1:2] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.50[,3:4] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.50) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.50[,3:4] <- NA #not estimated
qr.p.q.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.50) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.50[,1:2] <- t(sapply(qr.fit.50[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.50[,3:4] <- t(sapply(qr.fit.50[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.50[i,4]<- anova(qr.fit.50[["tau = 0.25"]][[i]],qr.fit.50[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.50 <- as.data.frame(rbind(two.step.p.50,ml.p.50,reml.p.50,qr.p.50))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.50 <- cbind(p.50,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.50 <- merge(coef.50, se.50, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.50 <- merge(est.50, p.50, by=c("repetition", "model"), sort=FALSE)
colnames(est.50)[colnames(est.50) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.50 <- reshape(est.50, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = estimands
)
# Include true values
est.long.50$true <- sapply(est.long.50$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.50 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.50 | two.step.gamma.fit.na.50),
                                  fit = two.step.fit.50[two.step.lm.fit.na.50 | two.step.gamma.fit.na.50]),
                 ML = cbind(repetition = which(ml.fit.na.50),
                            fit = ml.fit.50[ml.fit.na.50]),
                 REML = cbind(repetition = which(reml.fit.na.50),
                              fit = reml.fit.50[reml.fit.na.50]),
                 QR = cbind(repetition = which(qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50),
                            fit = qr.fit.50[qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50]))
saveRDS(error.50, file = "dataS5_50_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.50", "ml.fit.50", "reml.fit.50", "qr.fit.50"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.50$bias <- est.long.50$coef - est.long.50$true

# Coverage
# Construct CI first
conf.long.50 <- CI.fun(coef = est.long.50[,"coef"], se = est.long.50[,"se"], alpha=0.05)
est.long.50 <- cbind(est.long.50, conf.long.50)
# Check if true value is in the CI
est.long.50$cover <- (est.long.50$true >= est.long.50$ci_lo) & (est.long.50$true <= est.long.50$ci_hi)

# Rejection
reject <- est.long.50$p <= 0.05
est.long.50 <- cbind(est.long.50,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.50.mean <- aggregate(. ~ model+coef.type, data = est.long.50, na.action = na.pass,
                         FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.50, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.50.agg <- merge(est.50.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.50 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.50, na.action = na.pass, var, na.rm=T)
est.50.agg <- merge(est.50.agg, sampvar.50, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.50.agg$empse <- sqrt(est.50.agg$sampvar)

ms.50 <- multisimsum(
  data = est.long.50,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.50 <- summary(ms.50)$summ


# We need to calculate rejection MCSE on our own

est.50.agg$rejection_mcse <- sqrt(est.50.agg$reject* (1-est.50.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.50, file="S5_Nullinf_ms_summ_50.rds")
saveRDS(est.50.agg, file="S5_Nullinf_est_agg_50.rds")
saveRDS(est.long.50, file="S5_Nullinf_est_long_50.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

## 5.3 N-obs = 1000 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 1000
n_obs.1000 <- 1000

# Set nonrandom parameters
beta_0 <- 0.5
beta_1 <- 0.0
lambda_0 <- 0.5
lambda_1 <- 0
k <- 0.1 # 10% outlier
sigma <- 1

true.coef <- c(beta_0=beta_0, beta_1=beta_1, lambda_0=lambda_0, lambda_1=lambda_1)

#preallocate a list for the simulated data
simdata.1000 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.1000 <- .Random.seed

# Draw random data
simdata.1000 <- lapply(simdata.1000,function(r) {
  x1 <- runif(n_obs.1000*(1-k), 1, 10)
  normal <- rnorm(n_obs.1000*(1-k), mean = 0, sd = sigma)
  u <- sqrt(exp(lambda_0+lambda_1*x1))*normal
  y1 <- beta_0 + x1*beta_1 + u
  outlier <- mvrnorm(n_obs.1000*k, mu = c(8,beta_0+3.5+8*beta_1), Sigma = matrix(c(0.5,0,0,0.5),2,2))
  x1 <- c(x1,outlier[,1])
  y1 <- c(y1,outlier[,2])
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y1=y1,x1=x1,random=random)
})

#combine the randomnumberstates
random.1000 <- lapply(simdata.1000,function(data) data$random)
random.1000 <- c(list(random0.1000),random.1000)

## Determine true quantile regression values

tau <- c(0.25,0.5,0.75) #set quantiles to be estimated
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#partial deriviative of u: 1/2*lambda_0*sqrt(exp(lambda_0+lambda_1*x))

btau_0 <- beta_0 + sqrt(exp(lambda_0))*qnorm(tau, mean=0, sd=sigma)
btau_1 <- beta_1 + 1/2*lambda_1*sqrt(exp(lambda_0+lambda_1*mean(1:10)))*qnorm(tau, mean=0, sd=sigma)
true.btau <- rbind(btau_0,btau_1)
colnames(true.btau) <- c("tau = 0.25","tau = 0.5","tau = 0.75")

### Fit regressions #########################

# Fit two-step regression

two.step.fit.1000 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.1000 <- lapply(simdata.1000, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y1~data$x1),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x1, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.1000 <- !sapply(ml.fit.1000,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y1~data$x1,dformula = ~ data$x1+1, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.1000 <- !sapply(reml.fit.1000,is.model.fun,legalmod)

# Fit quantile regression

qr.fit.1000 <- matrix(data=NA,nrow=n_sim)
qr.fit.1000 <- lapply(tau,function(t) lapply(simdata.1000, function(data){
  tryCatch(rq(data$y1~data$x1, t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))
#check Fit
qr.t25.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.1000 <- !sapply(qr.fit.1000[["tau = 0.75"]],is.model.fun,legalmod)


#Build na.table
fit.na.1000 <- cbind(two.step.lm.fit.na.1000,two.step.gamma.fit.na.1000,ml.fit.na.1000,reml.fit.na.1000,qr.t25.fit.na.1000,qr.t50.fit.na.1000,qr.t75.fit.na.1000)
colSums(fit.na.1000)

#

# Delete simdata to save RAM
save(simdata.1000, file = "dataS5_1000_simdata.rds")
rm(list=c("simdata.1000"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.coef.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.coef.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.coef.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.coef.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.1000[,3:4] <- NA #not estimated
qr.coef.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.coef.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.coef.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
qr.coef.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Die Differenz zwischen beta_75 und beta_25 kann auch als Schätzer für lambda_1 betrachtet werden.
# Wegen IQR = 1,349*sigma, muss das Ergebnis aber skaliert werden:

qr.coef.1000[,3:4] <- (qr.coef.q.1000[,3:4]-qr.coef.q.1000[,1:2])/(qnorm(0.75,0,1)-qnorm(0.25,0,1))


# Combine estimates in Data.Frame
coef.1000 <- as.data.frame(rbind(two.step.coef.1000,ml.coef.1000,reml.coef.1000, qr.coef.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
coef.1000 <- cbind(coef.1000,repetition,model)

# Extract model SE

two.step.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.se.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.se.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.se.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,2]
  }else{
    c(NA,NA) #set coef to NA if model could not be estimated
  }
}))

qr.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.se.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.1000[,3:4] <- NA #not estimated
qr.se.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.se.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.se.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))
qr.se.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,2]
  }else{
    c(NA,NA) #set se to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.1000 <- as.data.frame(rbind(two.step.se.1000,ml.se.1000,reml.se.1000,qr.se.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
se.1000 <- cbind(se.1000,repetition,model)


# Extract p-values

two.step.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(two.step.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
two.step.p.1000[,1:2] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.1000[,3:4] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(ml.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
ml.p.1000[,1:2] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.1000[,3:4] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

reml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(reml.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
reml.p.1000[,1:2] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.1000[,3:4] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,4]
  }else{
    c(NA,NA) #set p-values to NA if model could not be estimated
  }
}))

qr.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.1000) <- c("beta_0","beta_1","lambda_0","lambda_1")
qr.p.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.5"]],function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.1000[,3:4] <- NA #not estimated
qr.p.q.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
colnames(qr.p.q.1000) <- c("beta_t25_0","beta_t25_x","beta_t75_0","beta_t75_x")
qr.p.q.1000[,1:2] <- t(sapply(qr.fit.1000[["tau = 0.25"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))
qr.p.q.1000[,3:4] <- t(sapply(qr.fit.1000[["tau = 0.75"]], function(fit){
  if(is.model.fun(fit,legalmod)){
    summary(fit, cov=TRUE)$coefficients[,4]
  }else{
    c(NA,NA) #pt p to NA if model could not be estimated
  }
}))

## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test überprüft ... und entspricht
## hier lambda_1=0
for (i in 1:n_sim){
  qr.p.1000[i,4]<- anova(qr.fit.1000[["tau = 0.25"]][[i]],qr.fit.1000[["tau = 0.75"]][[i]])$table$pvalue
}


# Build table of p-values of the coef
p.1000 <- as.data.frame(rbind(two.step.p.1000,ml.p.1000,reml.p.1000,qr.p.1000))
repetition <- rep(1:n_sim,4)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim),rep("QR",n_sim))
p.1000 <- cbind(p.1000,repetition,model)

# Combine all estimates
estimands <- c("beta_0","beta_1","lambda_0","lambda_1")

est.1000 <- merge(coef.1000, se.1000, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.1000 <- merge(est.1000, p.1000, by=c("repetition", "model"), sort=FALSE)
colnames(est.1000)[colnames(est.1000) %in% estimands] <- c("beta_0.p","beta_1.p","lambda_0.p","lambda_1.p")

est.long.1000 <- reshape(est.1000, idvar=c("repetition","model"), direction ="long",
                         varying = list(coef=c(3:6), se=c(7:10), p=c(11:14)),
                         v.names = c("coef","se","p"),
                         timevar = "coef.type",
                         times = estimands
)
# Include true values
est.long.1000$true <- sapply(est.long.1000$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.1000 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000),
                                    fit = two.step.fit.1000[two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000]),
                   ML = cbind(repetition = which(ml.fit.na.1000),
                              fit = ml.fit.1000[ml.fit.na.1000]),
                   REML = cbind(repetition = which(reml.fit.na.1000),
                                fit = reml.fit.1000[reml.fit.na.1000]),
                   QR = cbind(repetition = which(qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000),
                              fit = qr.fit.1000[qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000]))
saveRDS(error.1000, file = "dataS5_1000_errors.rds")


# Delete all fits to save RAM

rm(list=c("two.step.fit.1000", "ml.fit.1000", "reml.fit.1000", "qr.fit.1000"))

### Calculate Performance Measures and Monte-Carlo-Error #####


# Calculate some performance for single repetition first
# Bias calculation
est.long.1000$bias <- est.long.1000$coef - est.long.1000$true

# Coverage
# Construct CI first
conf.long.1000 <- CI.fun(coef = est.long.1000[,"coef"], se = est.long.1000[,"se"], alpha=0.05)
est.long.1000 <- cbind(est.long.1000, conf.long.1000)
# Check if true value is in the CI
est.long.1000$cover <- (est.long.1000$true >= est.long.1000$ci_lo) & (est.long.1000$true <= est.long.1000$ci_hi)

# Rejection
reject <- est.long.1000$p <= 0.05
est.long.1000 <- cbind(est.long.1000,reject)



# Average of bias, coverage and rejection is the estimate of the performance measure; na are omitted
est.1000.mean <- aggregate(. ~ model+coef.type, data = est.long.1000, na.action = na.pass,
                           FUN = mean, na.rm=T)


# Average model SE is the square root of the mean of the squared se
modelse <- aggregate(cbind(modelse=se) ~ model+coef.type, data = est.long.1000, na.action=na.pass,
                     FUN = function (x){
                       sqrt(mean(x^2, na.rm=T))
                     })
est.1000.agg <- merge(est.1000.mean, modelse, by=c("model","coef.type"), sort=FALSE)

# Calculate the sampling variance of all estimands
sampvar.1000 <- aggregate(cbind(sampvar=coef) ~ model+coef.type, data = est.long.1000, na.action = na.pass, var, na.rm=T)
est.1000.agg <- merge(est.1000.agg, sampvar.1000, by = c("model","coef.type"), sorte=FALSE)
# Calculate empirical se
est.1000.agg$empse <- sqrt(est.1000.agg$sampvar)

ms.1000 <- multisimsum(
  data = est.long.1000,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.1000 <- summary(ms.1000)$summ


# We need to calculate rejection MCSE on our own

est.1000.agg$rejection_mcse <- sqrt(est.1000.agg$reject* (1-est.1000.agg$reject) / n_sim)

### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.1000, file="S5_Nullinf_ms_summ_1000.rds")
saveRDS(est.1000.agg, file="S5_Nullinf_est_agg_1000.rds")
saveRDS(est.long.1000, file="S5_Nullinf_est_long_1000.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

# SZENARIO 6 Var-Exp --------------------------------------------------------------
# Multivariates, komplexes Mean-Modell und bivariates Scale-Modell
## 6.1 N-obs = 30 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 30
n_obs.30 <- 30

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

coef.name <- c("beta_0","beta_1","beta_2","beta_3","beta_4","lambda_0","lambda_1","lambda_2","lambda_3")
true.coef <- c(beta_0, beta_1, beta_2, beta_3, beta_4, lambda_0, lambda_1, lambda_2, lambda_3)
names(true.coef) <- coef.name

#preallocate a list for the simulated data
simdata.30 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.30 <- .Random.seed

# Draw random data
simdata.30 <- lapply(simdata.30,function(r) {
  x <- runif(n_obs.30, 1, 10)
  xx <- x^2
  z <- rbinom(n_obs.30,1,0.5)
  xz <- x*z
  normal <- rnorm(n_obs.30,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0 + lambda_1*x + lambda_2*z + lambda_3*x*z))*normal
  y <- beta_0+x*beta_1 + x^2*beta_2 + z*beta_3 + x*z*beta_4 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y=y, x=x, z=z, xx=xx, xz=xz, random=random)
})




#combine the randomnumberstates
random.30 <- lapply(simdata.30,function(data) data$random)
random.30 <- c(list(random0.30),random.30)

### Fit regressions #########################

# Fit two-step regression

two.step.fit.30 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.30 <- lapply(simdata.30, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y~data$x + data$xx + data$z + data$xz),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x + data$z + data$xz, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.30 <- !sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y~data$x + data$xx + data$z + data$xz,dformula = ~ data$x +  data$z + data$xz, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.30 <- !sapply(ml.fit.30,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.30 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.30 <- lapply(simdata.30, function(data){
  tryCatch(dglm(formula=data$y~data$x + data$xx + data$z + data$xz,dformula = ~ data$x +  data$z + data$xz, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.30 <- !sapply(reml.fit.30,is.model.fun,legalmod)


#Build na.table
fit.na.30 <- cbind(two.step.lm.fit.na.30,two.step.gamma.fit.na.30,ml.fit.na.30,reml.fit.na.30)
colSums(fit.na.30)

# 34 ML, 22 REML 52 Gamma

# Delete simdata to save RAM

save(simdata.30, file = "dataS6_30_simdata.rds")
rm(list=c("simdata.30"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(two.step.coef.30) <- coef.name
two.step.coef.30[,1:5] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.30[,6:9] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(ml.coef.30) <- coef.name
ml.coef.30[,1:5] <- t(sapply(ml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.30[,6:9] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.30 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(reml.coef.30) <- coef.name
reml.coef.30[,1:5] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.30[,6:9] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Combine estimates in Data.Frame
coef.30 <- as.data.frame(rbind(two.step.coef.30,ml.coef.30,reml.coef.30))
repetition <- rep(1:n_sim,3)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim))
coef.30 <- cbind(coef.30,repetition,model)

# Extract model SE

two.step.se.30 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(two.step.se.30) <- coef.name
two.step.se.30[,1:5] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.30[,6:9] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(ml.se.30) <- coef.name
ml.se.30[,1:5] <- t(sapply(1:length(ml.fit.30),function(fit){
  if(!any(is.na(ml.coef.30[fit,]))){
    as.numeric(coef(summary(ml.fit.30[[fit]]), dispersion=1)[,"Std. Error"])
  }else{
    as.numeric(c(NA,NA,NA,NA,NA)) #set coef to NA if model could not be estimated
  }
}))
ml.se.30[,6:9] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.30 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(reml.se.30) <- coef.name
reml.se.30[,1:5] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.30[,6:9] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.30 <- as.data.frame(rbind(two.step.se.30,ml.se.30,reml.se.30))
repetition <- rep(1:n_sim,3)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim))
se.30 <- cbind(se.30,repetition,model)


# Extract p-values

two.step.p.30 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(two.step.p.30) <- coef.name
two.step.p.30[,1:5] <- t(sapply(lapply(two.step.fit.30,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,"Pr(>|t|)"]
  }else{
    c(NA,NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.30[,6:9] <- t(sapply(lapply(two.step.fit.30,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(ml.p.30) <- coef.name
ml.p.30[,1:5] <- t(sapply(1:length(ml.fit.30),function(fit){
  if(!any(is.na(ml.coef.30[fit,]))){
    coef(summary(ml.fit.30[[fit]], dispersion=1))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.30[,6:9] <- t(sapply(lapply(ml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))


reml.p.30 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(reml.p.30) <- coef.name
reml.p.30[,1:5] <- t(sapply(reml.fit.30,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.30[,6:9] <- t(sapply(lapply(reml.fit.30,getElement,"dispersion.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Build table of p-values of the coef
p.30 <- as.data.frame(rbind(two.step.p.30,ml.p.30,reml.p.30))
repetition <- rep(1:n_sim,3)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim))
p.30 <- cbind(p.30,repetition,model)

# Combine all estimates

est.30 <- merge(coef.30, se.30, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.30 <- merge(est.30, p.30, by=c("repetition", "model"), sort=FALSE)
colnames(est.30)[colnames(est.30) %in% coef.name] <- c("beta_0.p","beta_1.p","beta_2.p","beta_3.p","beta_4.p","lambda_0.p","lambda_1.p","lambda_2.p","lambda_3.p")

est.long.30 <- reshape(est.30, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:11), se=c(12:20), p=c(21:29)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = coef.name
)
# Include true values
est.long.30$true <- sapply(est.long.30$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.30 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.30 | two.step.gamma.fit.na.30),
                                  fit = two.step.fit.30[two.step.lm.fit.na.30 | two.step.gamma.fit.na.30]),
                 ML = cbind(repetition = which(ml.fit.na.30),
                            fit = ml.fit.30[ml.fit.na.30]),
                 REML = cbind(repetition = which(reml.fit.na.30),
                              fit = reml.fit.30[reml.fit.na.30]))
saveRDS(error.30, file = "dataS6_30_errors.rds")


# Delete all fits to save RAM

rm(list=c("two.step.fit.30", "ml.fit.30", "reml.fit.30"))

### Calculate Performance Measures and Monte-Carlo-Error #####

ms.30 <- multisimsum(
  data = est.long.30,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.30 <- summary(ms.30)$summ


### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.30, file="S6_Exp_ms_summ_30.rds")
saveRDS(est.long.30, file="S6_Exp_est_long_30.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

## 6.2 N-obs = 50 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 50
n_obs.50 <- 50

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

coef.name <- c("beta_0","beta_1","beta_2","beta_3","beta_4","lambda_0","lambda_1","lambda_2","lambda_3")
true.coef <- c(beta_0, beta_1, beta_2, beta_3, beta_4, lambda_0, lambda_1, lambda_2, lambda_3)
names(true.coef) <- coef.name

#preallocate a list for the simulated data
simdata.50 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.50 <- .Random.seed

# Draw random data
simdata.50 <- lapply(simdata.50,function(r) {
  x <- runif(n_obs.50, 1, 10)
  xx <- x^2
  z <- rbinom(n_obs.50,1,0.5)
  xz <- x*z
  normal <- rnorm(n_obs.50,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0 + lambda_1*x + lambda_2*z + lambda_3*x*z))*normal
  y <- beta_0+x*beta_1 + x^2*beta_2 + z*beta_3 + x*z*beta_4 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y=y, x=x, z=z, xx=xx, xz=xz, random=random)
})




#combine the randomnumberstates
random.50 <- lapply(simdata.50,function(data) data$random) 
random.50 <- c(list(random0.50),random.50)

### Fit regressions #########################

# Fit two-step regression

two.step.fit.50 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.50 <- lapply(simdata.50, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y~data$x + data$xx + data$z + data$xz),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x + data$z + data$xz, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.50 <- !sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y~data$x + data$xx + data$z + data$xz,dformula = ~ data$x +  data$z + data$xz, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.50 <- !sapply(ml.fit.50,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.50 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.50 <- lapply(simdata.50, function(data){
  tryCatch(dglm(formula=data$y~data$x + data$xx + data$z + data$xz,dformula = ~ data$x +  data$z + data$xz, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.50 <- !sapply(reml.fit.50,is.model.fun,legalmod)


#Build na.table
fit.na.50 <- cbind(two.step.lm.fit.na.50,two.step.gamma.fit.na.50,ml.fit.na.50,reml.fit.na.50)
colSums(fit.na.50) 

# ML 17, REML 17, Gamma 38
# Delete simdata to save RAM

save(simdata.50, file = "dataS6_50_simdata.rds")
rm(list=c("simdata.50"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(two.step.coef.50) <- coef.name
two.step.coef.50[,1:5] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.50[,6:9] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(ml.coef.50) <- coef.name
ml.coef.50[,1:5] <- t(sapply(ml.fit.50,function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.50[,6:9] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.50 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(reml.coef.50) <- coef.name
reml.coef.50[,1:5] <- t(sapply(reml.fit.50,function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.50[,6:9] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Combine estimates in Data.Frame
coef.50 <- as.data.frame(rbind(two.step.coef.50,ml.coef.50,reml.coef.50))
repetition <- rep(1:n_sim,3)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim))
coef.50 <- cbind(coef.50,repetition,model)

# Extract model SE

two.step.se.50 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(two.step.se.50) <- coef.name
two.step.se.50[,1:5] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.50[,6:9] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(ml.se.50) <- coef.name
ml.se.50[,1:5] <- t(sapply(ml.fit.50,function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.50[,6:9] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.50 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(reml.se.50) <- coef.name
reml.se.50[,1:5] <- t(sapply(reml.fit.50,function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.50[,6:9] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.50 <- as.data.frame(rbind(two.step.se.50,ml.se.50,reml.se.50))
repetition <- rep(1:n_sim,3)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim))
se.50 <- cbind(se.50,repetition,model)


# Extract p-values

two.step.p.50 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(two.step.p.50) <- coef.name
two.step.p.50[,1:5] <- t(sapply(lapply(two.step.fit.50,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,"Pr(>|t|)"]
  }else{
    c(NA,NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.50[,6:9] <- t(sapply(lapply(two.step.fit.50,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(ml.p.50) <- coef.name
ml.p.50[,1:5] <- t(sapply(ml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.50[,6:9] <- t(sapply(lapply(ml.fit.50,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))


reml.p.50 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(reml.p.50) <- coef.name
reml.p.50[,1:5] <- t(sapply(reml.fit.50,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.50[,6:9] <- t(sapply(lapply(reml.fit.50,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Build table of p-values of the coef
p.50 <- as.data.frame(rbind(two.step.p.50,ml.p.50,reml.p.50))
repetition <- rep(1:n_sim,3)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim))
p.50 <- cbind(p.50,repetition,model)

# Combine all estimates

est.50 <- merge(coef.50, se.50, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.50 <- merge(est.50, p.50, by=c("repetition", "model"), sort=FALSE)
colnames(est.50)[colnames(est.50) %in% coef.name] <- c("beta_0.p","beta_1.p","beta_2.p","beta_3.p","beta_4.p","lambda_0.p","lambda_1.p","lambda_2.p","lambda_3.p")

est.long.50 <- reshape(est.50, idvar=c("repetition","model"), direction ="long",
                       varying = list(coef=c(3:11), se=c(12:20), p=c(21:29)),
                       v.names = c("coef","se","p"),
                       timevar = "coef.type",
                       times = coef.name
)
# Include true values
est.long.50$true <- sapply(est.long.50$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.50 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.50 | two.step.gamma.fit.na.50), 
                                  fit = two.step.fit.50[two.step.lm.fit.na.50 | two.step.gamma.fit.na.50]),
                 ML = cbind(repetition = which(ml.fit.na.50), 
                            fit = ml.fit.50[ml.fit.na.50]),
                 REML = cbind(repetition = which(reml.fit.na.50), 
                              fit = reml.fit.50[reml.fit.na.50]))
saveRDS(error.50, file = "dataS6_50_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.50", "ml.fit.50", "reml.fit.50"))

### Calculate Performance Measures and Monte-Carlo-Error #####

ms.50 <- multisimsum(
  data = est.long.50,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.50 <- summary(ms.50)$summ


### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.50, file="S6_Exp_ms_summ_50.rds")
saveRDS(est.long.50, file="S6_Exp_est_long_50.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

## 6.3 N-obs = 1000 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 1000
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

coef.name <- c("beta_0","beta_1","beta_2","beta_3","beta_4","lambda_0","lambda_1","lambda_2","lambda_3")
true.coef <- c(beta_0, beta_1, beta_2, beta_3, beta_4, lambda_0, lambda_1, lambda_2, lambda_3)
names(true.coef) <- coef.name

#preallocate a list for the simulated data
simdata.1000 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.1000 <- .Random.seed

# Draw random data
simdata.1000 <- lapply(simdata.1000,function(r) {
  x <- runif(n_obs.1000, 1, 10)
  xx <- x^2
  z <- rbinom(n_obs.1000,1,0.5)
  xz <- x*z
  normal <- rnorm(n_obs.1000,mean=0,sd=sigma)
  u <- sqrt(exp(lambda_0 + lambda_1*x + lambda_2*z + lambda_3*x*z))*normal
  y <- beta_0+x*beta_1 + x^2*beta_2 + z*beta_3 + x*z*beta_4 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y=y, x=x, z=z, xx=xx, xz=xz, random=random)
})




#combine the randomnumberstates
random.1000 <- lapply(simdata.1000,function(data) data$random) 
random.1000 <- c(list(random0.1000),random.1000)

### Fit regressions #########################

# Fit two-step regression

two.step.fit.1000 <- vector("list", length=n_sim) #preallocate to save calculation time
two.step.fit.1000 <- lapply(simdata.1000, function(data){
  lm.fit <-  tryCatch(lm(formula=data$y~data$x + data$xx + data$z + data$xz),
                      error=function(e) list(as.character(e))) #fit lm model of the mean and save error if impossible
  resid.lm <- residuals(lm.fit)^2
  gamma.fit <- tryCatch(glm(resid.lm~data$x + data$z + data$xz, family= Gamma (link=log)),
                        error=function(e) list(as.character(e))) #fit gamma glm model for the squared residuals and save error if impossible
  list(lm.fit = lm.fit,gamma.fit = gamma.fit)
})
#check Fit
two.step.lm.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),is.model.fun,legalmod)
two.step.gamma.fit.na.1000 <- !sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),is.model.fun,legalmod)

# Fit variance function regression with ml

ml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
ml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y~data$x + data$xx + data$z + data$xz,dformula = ~ data$x +  data$z + data$xz, family = gaussian(link = "identity"), dlink = "log", method="ml"),
           error=function(e) list(as.character(e))) #fit ml model for the variance function regression and save error if impossible
})
#check Fit
ml.fit.na.1000 <- !sapply(ml.fit.1000,is.model.fun,legalmod)

# Fit variance function regression with reml

reml.fit.1000 <- matrix(data=NA,nrow=n_sim) #preallocate to save calculation time
reml.fit.1000 <- lapply(simdata.1000, function(data){
  tryCatch(dglm(formula=data$y~data$x + data$xx + data$z + data$xz,dformula = ~ data$x +  data$z + data$xz, family = gaussian(link = "identity"), dlink = "log", method="reml"),
           error=function(e) list(as.character(e))) #fit reml model for the variance function regression and save error if impossible
})
#check Fit
reml.fit.na.1000 <- !sapply(reml.fit.1000,is.model.fun,legalmod)


#Build na.table
fit.na.1000 <- cbind(two.step.lm.fit.na.1000,two.step.gamma.fit.na.1000,ml.fit.na.1000,reml.fit.na.1000)
colSums(fit.na.1000) 

# Delete simdata to save RAM

save(simdata.1000, file = "dataS6_1000_simdata.rds")
rm(list=c("simdata.1000"))

### Extract model estimates ###########################

# Build coef.est

# Extract coeficients from the models if possible, otherwise set to NA
two.step.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(two.step.coef.1000) <- coef.name
two.step.coef.1000[,1:5] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
two.step.coef.1000[,6:9] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

ml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(ml.coef.1000) <- coef.name
ml.coef.1000[,1:5] <- t(sapply(ml.fit.1000,function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.coef.1000[,6:9] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.coef.1000 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(reml.coef.1000) <- coef.name
reml.coef.1000[,1:5] <- t(sapply(reml.fit.1000,function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.coef.1000[,6:9] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(fit)
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Combine estimates in Data.Frame
coef.1000 <- as.data.frame(rbind(two.step.coef.1000,ml.coef.1000,reml.coef.1000))
repetition <- rep(1:n_sim,3)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim))
coef.1000 <- cbind(coef.1000,repetition,model)

# Extract model SE

two.step.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(two.step.se.1000) <- coef.name
two.step.se.1000[,1:5] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA,NA) #set se to NA if model could not be estimated
  }
}))
two.step.se.1000[,6:9] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA) #set se to NA if model could not be estimated
  }
}))

ml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(ml.se.1000) <- coef.name
ml.se.1000[,1:5] <- t(sapply(ml.fit.1000,function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
ml.se.1000[,6:9] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

reml.se.1000 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(reml.se.1000) <- coef.name
reml.se.1000[,1:5] <- t(sapply(reml.fit.1000,function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit), dispersion=1)[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))
reml.se.1000[,6:9] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Std. Error"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

#build table of estimated se of the coef
se.1000 <- as.data.frame(rbind(two.step.se.1000,ml.se.1000,reml.se.1000))
repetition <- rep(1:n_sim,3)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim))
se.1000 <- cbind(se.1000,repetition,model)


# Extract p-values

two.step.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(two.step.p.1000) <- coef.name
two.step.p.1000[,1:5] <- t(sapply(lapply(two.step.fit.1000,getElement,"lm.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit))[,"Pr(>|t|)"]
  }else{
    c(NA,NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))
two.step.p.1000[,6:9] <- t(sapply(lapply(two.step.fit.1000,getElement,"gamma.fit"),function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))

ml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(ml.p.1000) <- coef.name
ml.p.1000[,1:5] <- t(sapply(ml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))
ml.p.1000[,6:9] <- t(sapply(lapply(ml.fit.1000,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))


reml.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=9)
colnames(reml.p.1000) <- coef.name
reml.p.1000[,1:5] <- t(sapply(reml.fit.1000,function(fit){
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=1))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA,NA) #set p-values to NA if model could not be estimated
  }
}))
reml.p.1000[,6:9] <- t(sapply(lapply(reml.fit.1000,getElement,"dispersion.fit"),function(fit){ 
  if(is.model.fun(fit,legalmod)){
    coef(summary(fit, dispersion=2))[,"Pr(>|z|)"]
  }else{
    c(NA,NA,NA,NA) #set coef to NA if model could not be estimated
  }
}))

# Build table of p-values of the coef
p.1000 <- as.data.frame(rbind(two.step.p.1000,ml.p.1000,reml.p.1000))
repetition <- rep(1:n_sim,3)
model <- c(rep("Two-Step",n_sim),rep("ML",n_sim),rep("REML",n_sim))
p.1000 <- cbind(p.1000,repetition,model)

# Combine all estimates

est.1000 <- merge(coef.1000, se.1000, by=c("repetition", "model"),suffixes=c(".coef",".se"), sort=FALSE)
est.1000 <- merge(est.1000, p.1000, by=c("repetition", "model"), sort=FALSE)
colnames(est.1000)[colnames(est.1000) %in% coef.name] <- c("beta_0.p","beta_1.p","beta_2.p","beta_3.p","beta_4.p","lambda_0.p","lambda_1.p","lambda_2.p","lambda_3.p")

est.long.1000 <- reshape(est.1000, idvar=c("repetition","model"), direction ="long",
                         varying = list(coef=c(3:11), se=c(12:20), p=c(21:29)),
                         v.names = c("coef","se","p"),
                         timevar = "coef.type",
                         times = coef.name
)
# Include true values
est.long.1000$true <- sapply(est.long.1000$coef.type, function(vec){
  true.coef[which(vec == names(true.coef))] #insert true-value depending on coef.type
}, USE.NAMES = FALSE
)

# Save error fits together with repetition number

error.1000 <- list(Two.Step = cbind(repetition = which(two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000), 
                                    fit = two.step.fit.1000[two.step.lm.fit.na.1000 | two.step.gamma.fit.na.1000]),
                   ML = cbind(repetition = which(ml.fit.na.1000), 
                              fit = ml.fit.1000[ml.fit.na.1000]),
                   REML = cbind(repetition = which(reml.fit.na.1000), 
                                fit = reml.fit.1000[reml.fit.na.1000]))
saveRDS(error.1000, file = "dataS6_1000_errors.rds")

# Delete all fits to save RAM

rm(list=c("two.step.fit.1000", "ml.fit.1000", "reml.fit.1000"))

### Calculate Performance Measures and Monte-Carlo-Error #####

ms.1000 <- multisimsum(
  data = est.long.1000,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "model",
  ref="Two-Step"
)

ms.summ.1000 <- summary(ms.1000)$summ


### Export data ###------------------------


# export data.frames for plotting and creating tables
saveRDS(ms.summ.1000, file="S6_Exp_ms_summ_1000.rds")
saveRDS(est.long.1000, file="S6_Exp_est_long_1000.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

# SZENARIO 7 Streuung Linear --------------------------------------------------------------
# Multivariate complex Mean und Scale data generation process
## 7.1 N-obs = 30 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 30
n_obs.30 <- 30

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

coef.name <- c("beta_0","beta_1","beta_2","beta_3","beta_4","lambda_0","lambda_1","lambda_2","lambda_3")
true.coef <- c(beta_0, beta_1, beta_2, beta_3, beta_4, lambda_0, lambda_1, lambda_2, lambda_3)
names(true.coef) <- coef.name

#preallocate a list for the simulated data
simdata.30 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.30 <- .Random.seed

# Draw random data
simdata.30 <- lapply(simdata.30,function(r) {
  x <- runif(n_obs.30, 1, 10)
  xx <- x^2
  z <- rbinom(n_obs.30,1,0.5)
  xz <- x*z
  normal <- rnorm(n_obs.30,mean=0,sd=sigma)
  u <- (lambda_0+lambda_1*x+lambda_2*z+lambda_3*x*z)*normal
  y <- beta_0+x*beta_1 + x^2*beta_2 + z*beta_3 + x*z*beta_4 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y=y, x=x, z=z, xx=xx, z=z, xz=xz, random=random)
})

#combine the randomnumberstates
random.30 <- lapply(simdata.30,function(data) data$random) 
random.30 <- c(list(random0.30),random.30)

# Derive true beta_tau values
tau <- c(0.25, 0.5, 0.75)
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#Median is equal to mean for normal distributed variables
true_betamed <- c(beta_0, beta_1, beta_2, beta_3, beta_4)
# for nomal distributed variables and variables that effect both mean and scale
# true beta_tau is equal to beta_median + lambda*z_tau, 
# where z is equal to the inverse of the normal cdf at the value of tau
# as x^2 is not part of the dispersion model, the true beta is also not affected (=0)
true_betatau <- true_betamed + t(qnorm(tau, 0, 1) %*% t(c(lambda_0, lambda_1, 0, lambda_2, lambda_3)))
true_betatau <- data.frame(true_betatau, row.names = c("(Intercept)", "x", "xx", "z", "xz"),
                           stringsAsFactors = FALSE)
names(true_betatau) <- c("Tau25", "Tau50", "Tau75")


### Fit regressions #########################

# Fit qr model directly

qr.fit.30 <- vector(mode = "list",length=n_sim)
qr.fit.30 <- lapply(simdata.30, function(d){
  
  tryCatch(rq(y~x + xx + z + xz, 
              data=data.frame(d[c("y","x","xx","z","xz")]),
              tau=tau),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
})

#check Fit
qr.fit.na <- !sapply(qr.fit.30, is.model.fun, legalmod)

# Fit qr model seperately for each tau

# qr.sep.fit.30 <- matrix(data=NA,nrow=n_sim)
qr.sep.fit.30 <- lapply(tau,function(t) lapply(simdata.30, function(d){
  tryCatch(rq(y~x + xx + z + xz, 
              data=data.frame(d[c("y","x","xx","z","xz")]), 
              tau=t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))

#check Fit
qr.t25.fit.na.30 <- !sapply(qr.sep.fit.30[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.30 <- !sapply(qr.sep.fit.30[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.30 <- !sapply(qr.sep.fit.30[["tau = 0.75"]],is.model.fun,legalmod)

save(simdata.30, file = "dataS7_30_simdata.rds")
rm(list=c("simdata.30"))

### Extract model estimates ###########################

# Build coef.est

qr.coef.30 <- lapply(1:length(qr.fit.30), function(r){
  repetition <- rep(r, nrow(true_betatau))
  coef.type <- rownames(coef(qr.fit.30[[r]]))
  if(is.model.fun(qr.fit.30[[r]], legalmod)){
    data.frame(coef.type, repetition, coef(qr.fit.30[[r]]))
  }else{
    data.frame(coef.type, repetition, matrix(NA, nrow= nrow(true_betatau), ncol= ncol(true_betatau))) #set coef to NA if model could not be estimated
  }
})

qr.coef.30 <- do.call(rbind, qr.coef.30)
names(qr.coef.30)[3:5] <- names(true_betatau)


# Extract model SE

qr.se.30 <- lapply(1:length(qr.fit.30), function(r){
  repetition <- rep(r, nrow(true_betatau))
  coef.type <- rownames(coef(qr.fit.30[[r]]))
  if(is.model.fun(qr.fit.30[[r]], legalmod)){
    se <- data.frame(sapply(lapply(summary(qr.fit.30[[r]], cov=TRUE), "[[", "coefficients"),function(m) m[,"Std. Error"]))
    names(se) <- names(true_betatau)
    data.frame(coef.type, repetition, se)
  }else{
    data.frame(coef.type, repetition, matrix(NA, nrow= nrow(true_betatau), ncol= ncol(true_betatau))) #set coef to NA if model could not be estimated
  }
})

qr.se.30 <- do.call(rbind, qr.se.30)


qr.diff.coef.30 <- qr.coef.30[,c("coef.type","repetition")]
qr.diff.coef.30$diff <- qr.coef.30[["Tau75"]]-qr.coef.30[["Tau25"]]


## Calculate p-values for difference in beta_tau
## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test getestet ...  

qr.diff.p.30 <- matrix(data=NA, nrow=n_sim,ncol=4)
for (i in 1:n_sim){
  qr.diff.p.30[i,] <- anova(qr.sep.fit.30[["tau = 0.25"]][[i]],
                            qr.sep.fit.30[["tau = 0.75"]][[i]],
                            joint="FALSE")$table$pvalue
}
qr.diff.p.30 <- data.frame(1:n_sim, qr.diff.p.30)
colnames(qr.diff.p.30) <- c("repetition",row.names(true_betatau)[-1])

qr.diff.p.long.30 <- reshape(qr.diff.p.30, direction ="long", 
                             idvar=c("repetition"),
                             varying= c("x","xx","z","xz"),
                             v.names = "p",
                             timevar = "coef.type",
                             times = c("x","xx","z","xz"))

qr.diff.long.30 <- merge(qr.diff.p.long.30, qr.diff.coef.30, by=c("coef.type", "repetition"), all.y=FALSE, sort=FALSE)
qr.diff.long.30$reject <-  qr.diff.long.30$p <= 0.05
true_diff <- true_betatau["Tau75"]-true_betatau["Tau25"]

# Include true values
qr.diff.long.30$true <- apply(qr.diff.long.30, 1, function(v){
  r <- match(v["coef.type"],row.names(true_diff))
  true_diff[r,]
}
)

qr.diff.agg.long.30 <-aggregate(.~ coef.type, data=qr.diff.long.30,
                                FUN = mean, na.rm=TRUE)
# Calculate rejection_mcse
qr.diff.agg.long.30$rejection_mcse <- sqrt(qr.diff.agg.long.30$reject* (1-qr.diff.agg.long.30$reject) / n_sim)


# Combine all estimates

qr.est.30 <- merge(qr.coef.30, qr.se.30, by=c("coef.type", "repetition"),suffixes=c(".coef",".se"), sort=FALSE)
qr.est.long.30 <- reshape(qr.est.30, direction ="long", 
                          idvar=c("coef.type","repetition"),
                          varying = list(coef=c(3:5), se=c(6:8)),
                          v.names = c("coef","se"),
                          timevar = "tau",
                          times = names(true_betatau))
# Include true values
qr.est.long.30$true <- apply(qr.est.long.30, 1, function(v){
  r <- match(v["coef.type"],row.names(true_betatau))
  c <- match(v["tau"],names(true_betatau))
  true_betatau[r,c]
}
)

# Calculate confidence intervals for the coefficients, assuming normal theory

qr.conf.long.30 <- CI.fun(coef = qr.est.long.30[,"coef"], se = qr.est.long.30[,"se"], alpha=0.05)
qr.est.long.30 <- cbind(qr.est.long.30, qr.conf.long.30)

#check Fit
qr.t25.fit.na.30 <- !sapply(qr.sep.fit.30[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.30 <- !sapply(qr.sep.fit.30[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.30 <- !sapply(qr.sep.fit.30[["tau = 0.75"]],is.model.fun,legalmod)

# Save error fits together with repetition number

error.30 <- list(QR_combined = cbind(repetition = which(qr.fit.na), 
                                     fit = qr.fit.30[qr.fit.na]),
                 QR_seperate = cbind(repetition = which(qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30), 
                                     fit = qr.sep.fit.30[qr.t25.fit.na.30 | qr.t50.fit.na.30 | qr.t75.fit.na.30]))
saveRDS(error.30, file = "dataS7_30_errors.rds")

# Delete all fits to save RAM

rm(list=c("qr.fit.30"))

### Calculate Performance Measures


### Calculate Performance Measures and Monte-Carlo-Error #####

ms.30 <- multisimsum(
  data = qr.est.long.30,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "tau",
  ref="Tau50"
)

ms.summ.30 <- summary(ms.30)$summ

ms.diff.30 <- multisimsum(
  data = qr.diff.long.30,
  par = "coef.type",
  true = "true",
  estvarname = "diff"
)


ms.diff.summ.30 <- summary(ms.diff.30)$summ


### Export data ###------------------------


# Export data.frames for plotting and creating tables
saveRDS(ms.summ.30, file="S7_Lin_ms_summ_30.rds")
saveRDS(qr.est.long.30, file="S7_Lin_est_long_30.rds")
saveRDS(ms.diff.summ.30, file="S7_Lin_diff_ms_summ_30.rds")
saveRDS(qr.diff.agg.long.30, file="S7_Lin_est_agg_long_30.rds")

# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])


## 7.2 N-obs = 50 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 50
n_obs.50 <- 50

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

coef.name <- c("beta_0","beta_1","beta_2","beta_3","beta_4","lambda_0","lambda_1","lambda_2","lambda_3")
true.coef <- c(beta_0, beta_1, beta_2, beta_3, beta_4, lambda_0, lambda_1, lambda_2, lambda_3)
names(true.coef) <- coef.name

#preallocate a list for the simulated data
simdata.50 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.50 <- .Random.seed

# Draw random data
simdata.50 <- lapply(simdata.50,function(r) {
  x <- runif(n_obs.50, 1, 10)
  xx <- x^2
  z <- rbinom(n_obs.50,1,0.5)
  xz <- x*z
  normal <- rnorm(n_obs.50,mean=0,sd=sigma)
  u <- (lambda_0+lambda_1*x+lambda_2*z+lambda_3*x*z)*normal
  y <- beta_0+x*beta_1 + x^2*beta_2 + z*beta_3 + x*z*beta_4 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y=y, x=x, z=z, xx=xx, z=z, xz=xz, random=random)
})

#combine the randomnumberstates
random.50 <- lapply(simdata.50,function(data) data$random) 
random.50 <- c(list(random0.50),random.50)

# Derive true beta_tau values
tau <- c(0.25, 0.5, 0.75)
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#Median is equal to mean for normal distributed variables
true_betamed <- c(beta_0, beta_1, beta_2, beta_3, beta_4)
# for nomal distributed variables and variables that effect both mean and scale
# true beta_tau is equal to beta_median + lambda*z_tau, 
# where z is equal to the inverse of the normal cdf at the value of tau
# as x^2 is not part of the dispersion model, the true beta is also not affected (=0)
true_betatau <- true_betamed + t(qnorm(tau, 0, 1) %*% t(c(lambda_0, lambda_1, 0, lambda_2, lambda_3)))
true_betatau <- data.frame(true_betatau, row.names = c("(Intercept)", "x", "xx", "z", "xz"),
                           stringsAsFactors = FALSE)
names(true_betatau) <- c("Tau25", "Tau50", "Tau75")


### Fit regressions #########################

# Fit qr model directly

qr.fit.50 <- vector(mode = "list",length=n_sim)
qr.fit.50 <- lapply(simdata.50, function(d){
  
  tryCatch(rq(y~x + xx + z + xz, 
              data=data.frame(d[c("y","x","xx","z","xz")]),
              tau=tau),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
})

#check Fit
qr.fit.na <- !sapply(qr.fit.50, is.model.fun, legalmod)

# Fit qr model seperately for each tau

# qr.sep.fit.50 <- matrix(data=NA,nrow=n_sim)
qr.sep.fit.50 <- lapply(tau,function(t) lapply(simdata.50, function(d){
  tryCatch(rq(y~x + xx + z + xz, 
              data=data.frame(d[c("y","x","xx","z","xz")]), 
              tau=t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))

#check Fit
qr.t25.fit.na.50 <- !sapply(qr.sep.fit.50[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.50 <- !sapply(qr.sep.fit.50[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.50 <- !sapply(qr.sep.fit.50[["tau = 0.75"]],is.model.fun,legalmod)

save(simdata.50, file = "dataS7_50_simdata.rds")
rm(list=c("simdata.50"))

### Extract model estimates ###########################

# Build coef.est

qr.coef.50 <- lapply(1:length(qr.fit.50), function(r){
  repetition <- rep(r, nrow(true_betatau))
  coef.type <- rownames(coef(qr.fit.50[[r]]))
  if(is.model.fun(qr.fit.50[[r]], legalmod)){
    data.frame(coef.type, repetition, coef(qr.fit.50[[r]]))
  }else{
    data.frame(coef.type, repetition, matrix(NA, nrow= nrow(true_betatau), ncol= ncol(true_betatau))) #set coef to NA if model could not be estimated
  }
})

qr.coef.50 <- do.call(rbind, qr.coef.50)
names(qr.coef.50)[3:5] <- names(true_betatau)


# Extract model SE

qr.se.50 <- lapply(1:length(qr.fit.50), function(r){
  repetition <- rep(r, nrow(true_betatau))
  coef.type <- rownames(coef(qr.fit.50[[r]]))
  if(is.model.fun(qr.fit.50[[r]], legalmod)){
    se <- data.frame(sapply(lapply(summary(qr.fit.50[[r]], cov=TRUE), "[[", "coefficients"),function(m) m[,"Std. Error"]))
    names(se) <- names(true_betatau)
    data.frame(coef.type, repetition, se)
  }else{
    data.frame(coef.type, repetition, matrix(NA, nrow= nrow(true_betatau), ncol= ncol(true_betatau))) #set coef to NA if model could not be estimated
  }
})

qr.se.50 <- do.call(rbind, qr.se.50)


qr.diff.coef.50 <- qr.coef.50[,c("coef.type","repetition")]
qr.diff.coef.50$diff <- qr.coef.50[["Tau75"]]-qr.coef.50[["Tau25"]]


## Calculate p-values for difference in beta_tau
## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test getestet ...  

qr.diff.p.50 <- matrix(data=NA, nrow=n_sim,ncol=4)
for (i in 1:n_sim){
  qr.diff.p.50[i,] <- anova(qr.sep.fit.50[["tau = 0.25"]][[i]],
                            qr.sep.fit.50[["tau = 0.75"]][[i]],
                            joint="FALSE")$table$pvalue
}
qr.diff.p.50 <- data.frame(1:n_sim, qr.diff.p.50)
colnames(qr.diff.p.50) <- c("repetition",row.names(true_betatau)[-1])

qr.diff.p.long.50 <- reshape(qr.diff.p.50, direction ="long", 
                             idvar=c("repetition"),
                             varying= c("x","xx","z","xz"),
                             v.names = "p",
                             timevar = "coef.type",
                             times = c("x","xx","z","xz"))

qr.diff.long.50 <- merge(qr.diff.p.long.50, qr.diff.coef.50, by=c("coef.type", "repetition"), all.y=FALSE, sort=FALSE)
qr.diff.long.50$reject <-  qr.diff.long.50$p <= 0.05
true_diff <- true_betatau["Tau75"]-true_betatau["Tau25"]

# Include true values
qr.diff.long.50$true <- apply(qr.diff.long.50, 1, function(v){
  r <- match(v["coef.type"],row.names(true_diff))
  true_diff[r,]
}
)

qr.diff.agg.long.50 <-aggregate(.~ coef.type, data=qr.diff.long.50,
                                FUN = mean, na.rm=TRUE)
# Calculate rejection_mcse
qr.diff.agg.long.50$rejection_mcse <- sqrt(qr.diff.agg.long.50$reject* (1-qr.diff.agg.long.50$reject) / n_sim)

# Combine all estimates

qr.est.50 <- merge(qr.coef.50, qr.se.50, by=c("coef.type", "repetition"),suffixes=c(".coef",".se"), sort=FALSE)
qr.est.long.50 <- reshape(qr.est.50, direction ="long", 
                          idvar=c("coef.type","repetition"),
                          varying = list(coef=c(3:5), se=c(6:8)),
                          v.names = c("coef","se"),
                          timevar = "tau",
                          times = names(true_betatau))
# Include true values
qr.est.long.50$true <- apply(qr.est.long.50, 1, function(v){
  r <- match(v["coef.type"],row.names(true_betatau))
  c <- match(v["tau"],names(true_betatau))
  true_betatau[r,c]
}
)

# Calculate confidence intervals for the coefficients, assuming normal theory

qr.conf.long.50 <- CI.fun(coef = qr.est.long.50[,"coef"], se = qr.est.long.50[,"se"], alpha=0.05)
qr.est.long.50 <- cbind(qr.est.long.50, qr.conf.long.50)

# Save errors

error.50 <- list(QR_combined = cbind(repetition = which(qr.fit.na), 
                                     fit = qr.fit.50[qr.fit.na]),
                 QR_seperate = cbind(repetition = which(qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50), 
                                     fit = qr.sep.fit.50[qr.t25.fit.na.50 | qr.t50.fit.na.50 | qr.t75.fit.na.50]))
saveRDS(error.50, file = "dataS7_50_errors.rds")

# Delete all fits to save RAM

rm(list=c("qr.fit.50"))

### Calculate Performance Measures


### Calculate Performance Measures and Monte-Carlo-Error #####

ms.50 <- multisimsum(
  data = qr.est.long.50,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "tau",
  ref="Tau50"
)

ms.summ.50 <- summary(ms.50)$summ

ms.diff.50 <- multisimsum(
  data = qr.diff.long.50,
  par = "coef.type",
  true = "true",
  estvarname = "diff"
)


ms.diff.summ.50 <- summary(ms.diff.50)$summ


### Export data ###------------------------


# Export data.frames for plotting and creating tables
saveRDS(ms.summ.50, file="S7_Lin_ms_summ_50.rds")
saveRDS(qr.est.long.50, file="S7_Lin_est_long_50.rds")
saveRDS(ms.diff.summ.50, file="S7_Lin_diff_ms_summ_50.rds")
saveRDS(qr.diff.agg.long.50, file="S7_Lin_est_agg_long_50.rds")

# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

## 7.3 N-obs = 1000 --------------------------------------------------------------------

### Simulate data -----------------------------------

#set number of replications in the simulation n_sim
n_sim <- 5000

# Set the number of observations n_obs = 1000
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

coef.name <- c("beta_0","beta_1","beta_2","beta_3","beta_4","lambda_0","lambda_1","lambda_2","lambda_3")
true.coef <- c(beta_0, beta_1, beta_2, beta_3, beta_4, lambda_0, lambda_1, lambda_2, lambda_3)
names(true.coef) <- coef.name

#preallocate a list for the simulated data
simdata.1000 <- vector("list", length=n_sim)

# store RandomNumberstate before the first repition
random0.1000 <- .Random.seed

# Draw random data
simdata.1000 <- lapply(simdata.1000,function(r) {
  x <- runif(n_obs.1000, 1, 10)
  xx <- x^2
  z <- rbinom(n_obs.1000,1,0.5)
  xz <- x*z
  normal <- rnorm(n_obs.1000,mean=0,sd=sigma)
  u <- (lambda_0+lambda_1*x+lambda_2*z+lambda_3*x*z)*normal
  y <- beta_0+x*beta_1 + x^2*beta_2 + z*beta_3 + x*z*beta_4 + u
  random <- .Random.seed #store RandomNumberstate after each repetition
  list(y=y, x=x, z=z, xx=xx, z=z, xz=xz, random=random)
})

#combine the randomnumberstates
random.1000 <- lapply(simdata.1000,function(data) data$random) 
random.1000 <- c(list(random0.1000),random.1000)

# Derive true beta_tau values
tau <- c(0.25, 0.5, 0.75)
names(tau) = c("tau = 0.25","tau = 0.5","tau = 0.75")

#Median is equal to mean for normal distributed variables
true_betamed <- c(beta_0, beta_1, beta_2, beta_3, beta_4)
# for nomal distributed variables and variables that effect both mean and scale
# true beta_tau is equal to beta_median + lambda*z_tau, 
# where z is equal to the inverse of the normal cdf at the value of tau
# as x^2 is not part of the dispersion model, the true beta is also not affected (=0)
true_betatau <- true_betamed + t(qnorm(tau, 0, 1) %*% t(c(lambda_0, lambda_1, 0, lambda_2, lambda_3)))
true_betatau <- data.frame(true_betatau, row.names = c("(Intercept)", "x", "xx", "z", "xz"),
                           stringsAsFactors = FALSE)
names(true_betatau) <- c("Tau25", "Tau50", "Tau75")


### Fit regressions #########################

# Fit qr model directly

qr.fit.1000 <- vector(mode = "list",length=n_sim)
qr.fit.1000 <- lapply(simdata.1000, function(d){
  
  tryCatch(rq(y~x + xx + z + xz, 
              data=data.frame(d[c("y","x","xx","z","xz")]),
              tau=tau),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
})

#check Fit
qr.fit.na <- !sapply(qr.fit.1000, is.model.fun, legalmod)

# Fit qr model seperately for each tau

# qr.sep.fit.1000 <- matrix(data=NA,nrow=n_sim)
qr.sep.fit.1000 <- lapply(tau,function(t) lapply(simdata.1000, function(d){
  tryCatch(rq(y~x + xx + z + xz, 
              data=data.frame(d[c("y","x","xx","z","xz")]), 
              tau=t),
           error=function(e) list(as.character(e))) #fit qr model and save error if impossible
}))

#check Fit
qr.t25.fit.na.1000 <- !sapply(qr.sep.fit.1000[["tau = 0.25"]],is.model.fun,legalmod)
qr.t50.fit.na.1000 <- !sapply(qr.sep.fit.1000[["tau = 0.5"]],is.model.fun,legalmod)
qr.t75.fit.na.1000 <- !sapply(qr.sep.fit.1000[["tau = 0.75"]],is.model.fun,legalmod)


save(simdata.1000, file = "dataS7_1000_simdata.rds")
rm(list=c("simdata.1000"))

### Extract model estimates ###########################

# Build coef.est

qr.coef.1000 <- lapply(1:length(qr.fit.1000), function(r){
  repetition <- rep(r, nrow(true_betatau))
  coef.type <- rownames(coef(qr.fit.1000[[r]]))
  if(is.model.fun(qr.fit.1000[[r]], legalmod)){
    data.frame(coef.type, repetition, coef(qr.fit.1000[[r]]))
  }else{
    data.frame(coef.type, repetition, matrix(NA, nrow= nrow(true_betatau), ncol= ncol(true_betatau))) #set coef to NA if model could not be estimated
  }
})

qr.coef.1000 <- do.call(rbind, qr.coef.1000)
names(qr.coef.1000)[3:5] <- names(true_betatau)


# Extract model SE

qr.se.1000 <- lapply(1:length(qr.fit.1000), function(r){
  repetition <- rep(r, nrow(true_betatau))
  coef.type <- rownames(coef(qr.fit.1000[[r]]))
  if(is.model.fun(qr.fit.1000[[r]], legalmod)){
    se <- data.frame(sapply(lapply(summary(qr.fit.1000[[r]], cov=TRUE), "[[", "coefficients"),function(m) m[,"Std. Error"]))
    names(se) <- names(true_betatau)
    data.frame(coef.type, repetition, se)
  }else{
    data.frame(coef.type, repetition, matrix(NA, nrow= nrow(true_betatau), ncol= ncol(true_betatau))) #set coef to NA if model could not be estimated
  }
})

qr.se.1000 <- do.call(rbind, qr.se.1000)


qr.diff.coef.1000 <- qr.coef.1000[,c("coef.type","repetition")]
qr.diff.coef.1000$diff <- qr.coef.1000[["Tau75"]]-qr.coef.1000[["Tau25"]]


## Calculate p-values for difference in beta_tau
## Die Nullhypothese beta_75=beta_25 wird über einen Anova-F-Test getestet ...  

qr.diff.p.1000 <- matrix(data=NA, nrow=n_sim,ncol=4)
for (i in 1:n_sim){
  qr.diff.p.1000[i,] <- anova(qr.sep.fit.1000[["tau = 0.25"]][[i]],
                              qr.sep.fit.1000[["tau = 0.75"]][[i]],
                              joint="FALSE")$table$pvalue
}
qr.diff.p.1000 <- data.frame(1:n_sim, qr.diff.p.1000)
colnames(qr.diff.p.1000) <- c("repetition",row.names(true_betatau)[-1])

qr.diff.p.long.1000 <- reshape(qr.diff.p.1000, direction ="long", 
                               idvar=c("repetition"),
                               varying= c("x","xx","z","xz"),
                               v.names = "p",
                               timevar = "coef.type",
                               times = c("x","xx","z","xz"))

qr.diff.long.1000 <- merge(qr.diff.p.long.1000, qr.diff.coef.1000, by=c("coef.type", "repetition"), all.y=FALSE, sort=FALSE)
qr.diff.long.1000$reject <-  qr.diff.long.1000$p <= 0.05

true_diff <- true_betatau["Tau75"]-true_betatau["Tau25"]

# Include true values
qr.diff.long.1000$true <- apply(qr.diff.long.1000, 1, function(v){
  r <- match(v["coef.type"],row.names(true_diff))
  true_diff[r,]
}
)

qr.diff.agg.long.1000 <-aggregate(.~ coef.type, data=qr.diff.long.1000,
                                  FUN = mean, na.rm=TRUE)
qr.diff.agg.long.1000$rejection_mcse <- sqrt(qr.diff.agg.long.1000$reject* (1-qr.diff.agg.long.1000$reject) / n_sim)

# Combine all estimates

qr.est.1000 <- merge(qr.coef.1000, qr.se.1000, by=c("coef.type", "repetition"),suffixes=c(".coef",".se"), sort=FALSE)
qr.est.long.1000 <- reshape(qr.est.1000, direction ="long", 
                            idvar=c("coef.type","repetition"),
                            varying = list(coef=c(3:5), se=c(6:8)),
                            v.names = c("coef","se"),
                            timevar = "tau",
                            times = names(true_betatau))
# Include true values
qr.est.long.1000$true <- apply(qr.est.long.1000, 1, function(v){
  r <- match(v["coef.type"],row.names(true_betatau))
  c <- match(v["tau"],names(true_betatau))
  true_betatau[r,c]
}
)

# Calculate confidence intervals for the coefficients, assuming normal theory

qr.conf.long.1000 <- CI.fun(coef = qr.est.long.1000[,"coef"], se = qr.est.long.1000[,"se"], alpha=0.05)
qr.est.long.1000 <- cbind(qr.est.long.1000, qr.conf.long.1000)

# Save error fits together with repetition number

error.1000 <- list(QR_combined = cbind(repetition = which(qr.fit.na), 
                                       fit = qr.fit.1000[qr.fit.na]),
                   QR_seperate = cbind(repetition = which(qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000), 
                                       fit = qr.sep.fit.1000[qr.t25.fit.na.1000 | qr.t50.fit.na.1000 | qr.t75.fit.na.1000]))
saveRDS(error.1000, file = "dataS7_1000_errors.rds")


# Delete all fits to save RAM

rm(list=c("qr.fit.1000"))

### Calculate Performance Measures


### Calculate Performance Measures and Monte-Carlo-Error #####

ms.1000 <- multisimsum(
  data = qr.est.long.1000,
  par = "coef.type",
  true = "true",
  estvarname = "coef",
  se = "se",
  methodvar = "tau",
  ref="Tau50"
)

ms.summ.1000 <- summary(ms.1000)$summ

ms.diff.1000 <- multisimsum(
  data = qr.diff.long.1000,
  par = "coef.type",
  true = "true",
  estvarname = "diff"
)


ms.diff.summ.1000 <- summary(ms.diff.1000)$summ


### Export data ###------------------------


# Export data.frames for plotting and creating tables
saveRDS(ms.summ.1000, file="S7_Lin_ms_summ_1000.rds")
saveRDS(qr.est.long.1000, file="S7_Lin_est_long_1000.rds")
saveRDS(ms.diff.summ.1000, file="S7_Lin_diff_ms_summ_1000.rds")
saveRDS(qr.diff.agg.long.1000, file="S7_Lin_est_agg_long_1000.rds")
# Delete everything except predefined variables and functions
rm(list=ls()[! (ls() %in% lsf.str() | ls() %in% c("legalmod"))])

