# before you run this code, run "read_data.R" script to load data
library(tidyverse)

N <- 10000

n <- 4000
p <- 800
rho <- 0.5
Sigma  <- toeplitz(rho^(0:(p-1)))
gamma2 <- 5
kappa <- 0.2

sigma_s = 4.7436
alpha_s = 1.4994
lambda_s = 3.0269

#-- true parameters
jnull <- scan(paste0(dir,'/section-5/param/jnull.txt'))
non_null <- scan(paste0(dir,'/section-5/param/non_null.txt'))


#-- obtain beta
beta <- numeric(p)
#find the beta value to get the desired gamma2
f <- function(x){
  beta[non_null] <- x
  
  return(t(beta)%*%(Sigma%*%beta)/n - gamma2)
}
beta_solve <- uniroot(f, c(0,300), check.conv=TRUE,extendInt = "upX")$root
beta[non_null] <- beta_solve

hat_alpha = param[,1];hat_sigma = param[,2];hat_lambda=param[,3]

# true tau
tau <- 1/sqrt(diag(solve(Sigma)))

#-- calculate tilde_beta using the three methods
tilde_beta_1 <- matrix(0, N, p)
tilde_beta_2 <- matrix(0, N, p)
tilde_beta_t <- matrix(0, N, p)

for(i in 1:N){
  tilde_beta_t[i,] <- tau*(beta_est[i,]-alpha_s*beta)/sigma_s # true parameter
  tilde_beta_1[i,] <- tau_hat_1[i,]*(beta_est[i,]-hat_alpha[i]*beta)/hat_sigma[i]
  tilde_beta_2[i,] <- tau_hat_2[i,]*(beta_est[i,]-hat_alpha[i]*beta)/hat_sigma[i]
}

#Table 3 : t-test
tilde_beta_1 <- beta_est/std_est # classical z-statistics
phi = 2-2*pnorm(abs(tilde_beta_t[,jnull]), 0, 1, lower.tail=TRUE)

alpha <- c(0.1,0.05,0.01,0.005)
output <- NULL
for(a in alpha){
  result <- c(a*100,mean(phi<a)*100, sqrt(mean(phi<a)*(1-mean(phi<a))/N)*100)
  
  output <- rbind(output, result)
}
output

# Table 4: marginal coverage of one coordinate
coverage_proportion <- tibble(
  prop = c(0.995, 0.99, 0.95, 0.90),
  normal_quantile = qnorm(0.5 + prop/2)
) %>% 
  mutate(
    cov_prop = normal_quantile %>% map_dbl(~ mean(tilde_beta_2[,jnull] < . & tilde_beta_2[,jnull] > -.)),
    sd = normal_quantile %>% map_dbl(~ sd(tilde_beta_2[,jnull] < . & tilde_beta_2[,jnull] > -.) / sqrt(10000))
)

coverage_proportion * 100

#plot an empirical distribution
qqnorm(tilde_beta_t[,jnull])
qqnorm(tilde_beta_1[,jnull])

qqline(tilde_beta_t[,jnull])
qqline(tilde_beta_1[,jnull])

# Table 5 : coverage proportion
tilde_beta <- tilde_beta_t

mean(rowMeans(abs(tilde_beta)<2.32))*100; sd(rowMeans(abs(tilde_beta)<2.32))/sqrt(N)*100 #98%
mean(rowMeans(abs(tilde_beta)<1.96))*100; sd(rowMeans(abs(tilde_beta)<1.96))/sqrt(N)*100 #95%
mean(rowMeans(abs(tilde_beta)<1.64))*100; sd(rowMeans(abs(tilde_beta)<1.64))/sqrt(N)*100 #90%
mean(rowMeans(abs(tilde_beta)<1.28))*100; sd(rowMeans(abs(tilde_beta)<1.28))/sqrt(N)*100 #80%

# Table 6: LLR
hist(pchisq(dev, lower.tail=FALSE, df=1))

adj_t <- kappa*sigma_s^2/lambda_s
adj_e <- kappa*hat_sigma[1:N]^2/hat_lambda[1:N]
adj_c <- 1
p.val=data.frame(pchisq(dev[1:N]/adj_e, lower.tail=FALSE, df=1))
phi <- p.val









