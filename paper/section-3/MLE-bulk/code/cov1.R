# This code is for random covariance case, the others are the same and thus not included

# job submission arguments
args = commandArgs(trailingOnly=TRUE)
ind <- args[1]

# -- Set parameters here
n <- 4000
kappa <- 0.2
p_sig <- 0.5
gamma2 <- 5

B <- 1000
# --

p <- n * kappa

# calculated paramters
alpha_s = 1.4994
sigma_s = 4.7436
lambda_s = 3.0269

# Read the covariance matrix that has been sampled before
Sigma <- as.matrix(read.table('/section-3/MLE-bulk/param/random/Sigma.txt', header=FALSE))
R <- chol(Sigma)
inv_Sigma <- 1/sqrt(diag(solve(Sigma)))
#sample beta
beta <- numeric(p)
non_null <- scan('/section-3/MLE-bulk/param/random/non_null.txt')
j_null <- scan('/section-3/MLE-bulk/param/random/jnull.txt')
beta_solve <- uniroot(function(x){ beta[non_null] <- x; return(t(beta)%*%(Sigma%*%beta)/n - gamma2);}, c(0,300), check.conv=TRUE,extendInt = "upX")$root
beta[non_null] <- beta_solve
print(abs( t(beta)%*%(Sigma%*%beta)/n - gamma2)) #check if desired signal strength
print(abs(t(beta)%*%beta)/n)

#confidence levels i will check
ci_level = c(0.99,0.98,0.95,0.90,0.80)
ci_thresh = c(2.676, 2.326, 1.96, 1.645, 1.281)
n_level = length(ci_level)

pc_t <- matrix(0, B, n_level)

fname = 'section-3/MLE-bulk/output/random/'

for(b in 1:B){
  X <- matrix(rnorm(n*p, 0, 1), n, p)%*%R/ sqrt(n); 
  eta <- X %*% beta; 
  mu <- exp(eta) / (1+exp(eta))
  y <- rbinom(n, 1, mu)
  
  cov_est <- 1/sqrt(diag(solve(t(X) %*% X)))/(sqrt(1-kappa)); 
  cov_est1[b] = cov_est[j_null]
  #fit a logistic regression
  fit <- glm(y~X+0, family=binomial)

  beta_est <- fit$coefficients; 
  std_store <- summary(fit)$coef[,2];  

  fit0 <- glm(y~X[,-j_null]+0, family=binomial)
  
  for(levels in 1:n_level){
    ci_t <-  cbind( fit$coefficients - ci_thresh[levels] * sigma_s/(inv_Sigma),fit$coefficients +  ci_thresh[levels] * sigma_s/(inv_Sigma))/alpha_s
    #compute coverage proportion
    pc_t[b, levels] <- mean((ci_t[,2] >= beta) & (ci_t[,1] <= beta))
  }
}

write.table(pc_t, paste0(fname,'pct',ind,'.txt'),col.names=FALSE, row.names=FALSE, append=FALSE)
