# Section 3.2.2 Condition on the regression coef
# compute the standard deviation of sqrt(n) (alpha_n - alpha_star) 

#-- job-submission arguments (please ignore)
args = commandArgs(trailingOnly=TRUE) 
kappa <- as.numeric(args[1])
gamma2 <- as.numeric(args[2])
n <- as.numeric(args[3])
ind <- as.numeric(args[4])
##-- 

#--- load beta
fileloc <- "/param/"
filename <- paste0(paste(n, kappa, gamma2, sep = "_"), ".txt")
beta <- scan(paste0(fileloc, filename))
p <- n * kappa

# covariance matrix
Sigma <- toeplitz(0.5^(0:(p-1)))
R <- chol(Sigma)

B <- 200 
alpha <- numeric(B)
for(b in 1:B){
  X <- matrix(rnorm(n*p, 0, 1), n, p) %*% R / sqrt(p); eta <- X %*% beta; mu <- 1 / (1+exp(-eta))
  y <- rbinom(n, 1, mu)
  
  #fit a logistic regression
  fit <- glm(y ~ X + 0, family=binomial)
  beta_hat <- fit$coefficients; 
  alpha[b] <- (t(beta_hat) %*% (Sigma %*% beta))[1,1] / p / gamma2 # compute alpah(n)
  
  cat(b, ",")
  if(b %% 20 == 0) cat("\n")
}

# store the computed alpha(n)
fileloc <- "/output/"
folder_name <- paste(kappa,gamma2,n,sep = "-")
write.table(alpha, paste0(fileloc, folder_name, "/alpha_", ind,".txt"),
            append = F, row.names = F, col.names = F)



