# Section 3.2.2 Condition on the regression coef 
# generate true coef beta when n varies 
# and sample non-null variables
library(tidyverse)
params <- data.frame(
  ratio = rep(c(0.15, 0.3, 0.5), each = 4),
  n = rep(c(1000, 2000, 4000, 8000), 3)
)

fileloc <- "/condition/param/"
kappa = 0.2
g = 2 # gamma 
for(i in 1:12){
  r = params[i,1]
  n <- params[i,2]
  p <- n * kappa
  Sigma <- toeplitz(0.5^(0:(p-1))) 
  tau <- 1/sqrt(diag(solve(Sigma)))
  
  # randomlu sample candidate non-nulls
  nonnull <- sample(1:p, 30, replace = F) 
  beta <- vector(length = p)
  diff <- numeric(30)
  # add non-nulls one by one untill we exceed the signal strength
  for(j in 1:30){
    beta[nonnull[1:j]] <- r * g / (tau[nonnull[1:j]] / sqrt(p))
    diff[j] <- (t(beta) %*% (Sigma %*% beta))[1,1] / p  - g^2
  }
  jj <- max(which(diff < 0))
  beta <- vector(length = p)
  # pick the last non-null coef value to achieve target gamma = 2
  beta[nonnull[1:(jj)]] <- r * g / (tau[nonnull[1:(jj)]]/ sqrt(p))
  sol <- uniroot(f = function(t) {beta[nonnull[jj+1]] <- t; (t(beta) %*% (Sigma %*% beta))[1,1] / p - g^2}, c(0, 10),
                 extendInt = "yes")
  beta[nonnull[jj+1]] <- sol$root
  
  cat("___kappa = ", kappa, ", ratio = ", r, "n = ",n,"\n")
  cat("Pick the non-null to be", nonnull[1], "\n")
  cat(beta[nonnull[1]] * tau[nonnull[1]] / g / sqrt(p), ",")
  cat(beta %*% (Sigma %*% beta) / p, "\n")
  
  filename <- paste(n, r, sep = "_")
  #    write.table(beta, file = paste0(fileloc, filename,".txt"),
  #                append = F, col.names = F, row.names = F)
}
