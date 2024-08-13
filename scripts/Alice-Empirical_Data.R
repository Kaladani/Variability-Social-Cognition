##############HEADER##############################
## author:    Daniel Gotthardt
## contact:   daniel.gotthardt@uni-hamburg.de / daniel.gotthardt@gmx.de
## file name: Alice-Empirical_data.R
## Context:   Theorizing and Testing Variability in Social Cognition
## Input:     -
## Output:    Simulated Dataset with Stereotype Valence and Prosocialness
## Summary:   This R-File generates empirical data for exposition in tutorial

################BODY##############################


# Clean up workspace -------------------------------------
rm(list=ls(all.names = TRUE))
gc()

# Helper functions ----------------------------------
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


# Simulate observed data ----

set.seed(13)

n.obs <- 250

X <- rnormt(n.obs, c(0,1), 0.5, 0.2)


X.mean <- rnormt(n.obs, c(0,1), 0.5, 0.2)

Z <- runif(n.obs, 0, 1)

# while we manipulated Z in the experimental design, we don't assume 
# that exemplar theory is true in this simulation, but that
# people have randomly varying estimation of the true mean

X.valence <- sapply(1:n.obs, function(i) {
  rnormc(1, c(0,1), X.mean[i], 0.001)
}
)

# X.valence <- X.mean
# we simulate the starting values as a normal distribution
# around a prosocialness of 0.4 with sd 0.2

pre <- rnormc(n.obs, c(0,1), 0.4, 0.2)

# defining the attractor 
Y.attractor <- 0.5 * X.valence + 0.4

# X.valence is also defining the strength of the attraction

Y.true <- pre + X.mean * (Y.attractor - pre)

# In the end some random measurement noise is added

Y.obs <- rnormc(n.obs, c(0,1), Y.true, 0.01)

# Everything is combined in a data.frame
mydata <- data.frame(X.mean = X.mean, 
                     X.valence = X.valence, 
                     Z = Z, 
                     Y.obs = Y.obs, 
                     Y.true = Y.true)

saveRDS(mydata, "data/Simulated_Empirical_Data.Rds")