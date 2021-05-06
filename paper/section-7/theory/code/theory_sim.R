# Section 7.1.1 Finite sample accuracy
location <- "/scratch/users/qzhao1/logistic/sim/cov/21apr9lrt/known/" 

#-- job-submission arguments (please ignore)
args = commandArgs(trailingOnly=TRUE) 
type <- args[1]
ind <- as.numeric(args[2])
#-- 

#--- Simulation parameters
n <- 4000
kappa <- 0.2
p <- 800 
gamma2 <- 5

B <- 1000 # Number of repetitions 

# Read the covariance matrix that has been sampled before
Sigma <- as.matrix(read.table(paste0(location, "param/Sigma_",type,".txt"), header=FALSE))
R <- chol(Sigma)
# read parameters
beta <- scan(paste0(location, "param/beta_",type,".txt"))
beta0 <- 1
# jnull = 478 # null coordinate to compute the LRT

#store 
mle <- matrix(0, B, p+1)
lrt <- numeric(B)
for(b in 1:B){
  X <- matrix(rnorm(n*p, 0, 1), n, p) %*% R / sqrt(p); eta <- X %*% beta + beta0; mu <- exp(eta) / (1+exp(eta))
  y <- rbinom(n, 1, mu)
  
  #fit a logistic regression
  fit <- glm(y ~ X, family=binomial)
  mle[b,] <- fit$coefficients;
  
  # the LRT is computed only for Sigma is ar(1) with rho = 0.5
  fit_small <- glm(y ~ X[ ,-jnull], family = binomial)
  lrt[b] <- anova(fit_small, fit, test = "LRT")$Deviance[2] 
}

# write.table(lrt, paste0(location,'output/', type, '/lrt_',ind,'.txt'), 
#            col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(mle, paste0(location,'output/', type, '/mle_',ind,'.txt'), 
           col.names=FALSE, row.names=FALSE, append=FALSE)
