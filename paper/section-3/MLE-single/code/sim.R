# Is this the file location
location <- "/section-3/"

#-- job-submission arguments (please ignore)
args = commandArgs(trailingOnly=TRUE) 
type <- args[1]
ind <- as.numeric(args[2])
##-- 

#--- Simulation parameters
n <- 4000
kappa <- 0.2
p <- 800 
p_sig <- 0.5
gamma2 <- 5

B <- 100 # Number of repetitions 

print(ind) # Debug

# Calculated solutions to the system of equations
alpha_s = 1.4994
sigma_s = 4.7436
lambda_s = 3.0269

# Read the covariance matrix that has been sampled before
Sigma <- as.matrix(read.table(paste0(location, "param/", type, "/Sigma.txt", header=FALSE)))
R <- chol(Sigma)
inv_Sigma <- 1/sqrt(diag(solve(Sigma)))
# Beta
beta <- numeric(p)
non_null <- scan(paste0(location, "param/", type, "/non_null.txt"))
j_nonnull <- scan(paste0(location, "param/", type, "/j_nonnull.txt"))
beta_solve <- uniroot(function(x){ beta[non_null] <- x; return(t(beta)%*%(Sigma%*%beta)/n - gamma2);}, c(0,300), check.conv=TRUE,extendInt = "upX")$root
beta[non_null] <- beta_solve
print(abs( t(beta)%*%(Sigma%*%beta)/n - gamma2)) # Check if desired signal strength
print(abs(t(beta)%*%beta)/n)

#store 
beta_est <- matrix(0, B, 10)

for(b in 1:B){
  X <- matrix(rnorm(n*p, 0, 1), n, p) %*% R / sqrt(n); eta <- X %*% beta; mu <- exp(eta) / (1+exp(eta))
  y <- rbinom(n, 1, mu)
  
  #fit a logistic regression
  fit <- glm(y~X+0, family=binomial)
  
  beta_est[b,] <- fit$coefficients[j_nonnull]; 
}

write.table(beta_est, paste0(location,'/output/', type, '/beta_est_',ind,'.txt'), 
            col.names=FALSE, row.names=FALSE, append=FALSE)


