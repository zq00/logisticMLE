# 7.2. Estimating model parameters
# Compute MLE
library(glmhd)

location <- "/unknown/" 

#-- job-submission arguments (please ignore)
args = commandArgs(trailingOnly=TRUE) 
type <- args[1]
ind <- as.numeric(args[2])
##-- 

#--- Simulation parameters
n <- 4000
kappa <- 0.2
p <- 800 

B <- 15 # Number of repetitions 

# Read the covariance matrix that has been sampled before
Sigma <- as.matrix(read.table(paste0(location, "param/Sigma_",type,".txt"), header=FALSE))
R <- chol(Sigma)
# read parameters
beta <- scan(paste0(location, "param/beta_",type,".txt"))
beta0 <- 1
# jnull = 478 # a null coef. for to compute the LRT

# mle 
mle <- matrix(0, B, p+1)
# estimated intercept and signal strength
estimatedGamma <- matrix(0, B, 2)
# estimated std
sigmahat <- matrix(0, B, p)
# estimated param
estimatedparam <- matrix(0, B, 4)
# lrt
lrt <- numeric(B)

for(b in 1:B){
  cat("--",b,"\n")
  X <- matrix(rnorm(n*p, 0, 1), n, p) %*% R / sqrt(p); eta <- X %*% beta + beta0; mu <- exp(eta) / (1+exp(eta))
  y <- rbinom(n, 1, mu)
  
  #fit a logistic regression
  fit <- glm(y ~ X, family=binomial, x = T, y = T)
  fit_adj <- adjust_glm(fit, verbose = T, echo = T)
  
  mle[b,] <- fit$coefficients; 
  estimatedGamma[b, ] <- c(fit_adj$intercept, fit_adj$gamma_hat)
  sigmahat[b, ] <- fit_adj$std_adj
  estimatedparam[b, ] <- fit_adj$param
  
  # compute the LRT only for covariance is AR(1) with rho = 0.5
  fit_small <- glm(y ~ X[ ,-jnull], family = binomial)
  lrt[b] <- anova(fit_small, fit, test = "LRT")$Deviance[2]
}

write.table(mle, paste0(location,'output/', type, '/mle_',ind,'.txt'), 
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(estimatedGamma, paste0(location,'output/', type, '/gammahat_',ind,'.txt'), 
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(sigmahat, paste0(location,'output/', type, '/sigmahat_',ind,'.txt'), 
            col.names=FALSE, row.names=FALSE, append=FALSE)

write.table(lrt, paste0(location,'/output/', type, '/lrt_',ind,'.txt'), 
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(estimatedparam, paste0(location,'/output/', type, '/param_',ind,'.txt'), 
            col.names=FALSE, row.names=FALSE, append=FALSE)
