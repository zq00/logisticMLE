# Section 3.2.2 Condition on the regression coef
# Figure 3 Simulation code 
# Compute the MLE for a non-null variable

#-- job-submission arguments (please ignore)
args = commandArgs(trailingOnly=TRUE) 
k <- as.numeric(args[1])
ind <- as.numeric(args[2])
##-- 

#--- Simulation parameters
# ratio - tau_j beta_j / gamma for j the non-null coord
# n - number of observations
# nonnull - the index of the non-null coord
params <- data.frame(
  ratio = rep(c(0.15, 0.3, 0.5), each = 4),
  n = rep(c(1000, 2000, 4000, 8000), 3),
  nonnull = c(80, 106, 133, 93,
              195, 30, 592, 825,
              3, 162, 793, 1576)
)

cat(k, "\n")
kappa <-  0.2
gamma <-  2
r <- params[k, 1]
n <- params[k, 2]
nonnull <- params[k, 3]
p = n * kappa

# read true beta
betaloc <- "/scratch/users/qzhao1/logistic/sim/cov/21apr27cond4/param/"
filename <- paste0(paste(n, r, sep = "_"), ".txt")
beta <- scan(paste0(betaloc, filename))

# covariance matrix
Sigma <- toeplitz(0.5^(0:(p-1)))
R <- chol(Sigma)

# check the signal strength and tau_j beta_j / gamma is as we want
cat( (t(beta) %*% (Sigma %*% beta))[1,1] / p)
cat(beta[nonnull] * 0.7745967 / gamma /sqrt(p)) # tau_j = 0.7745, divide by sqrt(p) because we scale X to have variance 1/p for each variable

B <- 1000  # repeat 1000 times. we then repeat each simulation 10 times. 
mle <- numeric(B)
for(b in 1:B){
  X <- matrix(rnorm(n*p, 0, 1), n, p) %*% R / sqrt(p); eta <- X %*% beta; mu <- 1 / (1+exp(-eta))
  y <- rbinom(n, 1, mu)
  
  #fit a logistic regression
  fit <- glm(y ~ X + 0, family=binomial)
  beta_hat <- fit$coefficients; 
  mle[b] <- beta_hat[nonnull]
  
  cat(b, ",")
  if(b %% 20 == 0) cat("\n")
}

fileloc <- "/condition/output/" # output location
folder_name <- paste(n, r, sep = "-")
write.table(mle, paste0(fileloc, folder_name, "/mle_", ind,".txt"),
            append = F, row.names = F, col.names = F)



