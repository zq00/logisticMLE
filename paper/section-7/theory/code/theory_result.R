# Section 7.1.1. Finite sample accuracy
# Results in Table 10 (I), 11 (I) and supplementary material Table 1 - 3
library(tidyverse)
library(glmhd)

type = "identity" # input the covariance type
fileloc <- "/theory/"
filenames <- list.files(paste0(fileloc, "output/", type , "/"))

# read the mle
beta_hat <- filenames %>% 
  map(~ read.table(paste0(fileloc,"output/", type , "/", .))) %>% 
  bind_rows()

N <- 100000
beta_hat <- beta_hat[1:N, ]

# Parameters
n <- 4000
p <- 800
gamma2 <- 5

# Covariance and true coef.
Sigma <- as.matrix(read.table(
  paste0(fileloc, "param/Sigma_", type,".txt")))

beta <- scan(paste0(fileloc, "param/beta_", type,".txt"))
beta0 <- 1
tau <- sqrt(1/diag(solve(Sigma))) # tau_j

# theoretical parameters
# params <- find_param(kappa = 0.2, gamma = sqrt(5), beta0 = 1)
alpha_s = 1.558095
sigma_s = 2.305896

# 1. Coverage proportion of a single non-null
# Table 10 Column I and supplementary table 1 
j <- sample(which(beta != 0), 1) # randomly pick a non-null
beta_adj <- (beta_hat[,j+1] - alpha_s * beta[j]) * tau[j] / sigma_s # adjusted beta_hat

# normal quantile plot 
data <- tibble(
  index = 1:N,
  beta_hat = sort(beta_adj)
) %>% 
  mutate(
    quantile = (index - 0.5) / N,
    normal_quantile = qnorm(quantile, 0, 1)
  ) 

# Coverage proportion
coverage_proportion <- tibble(
  prop = c(0.99, 0.98, 0.95, 0.90, 0.8),
  normal_quantile = qnorm(0.5 + prop/2)
) %>% 
  mutate(
    cov_prop = normal_quantile %>% map_dbl(~ mean(beta_adj < . & beta_adj > -.)),
    sd = normal_quantile %>% map_dbl(~ sd(beta_adj < . & beta_adj > -.) / sqrt(N))
  ) 

coverage_proportion*100

# 2. Coverage proportion of all the variables (bulk)
# Table 11 Column I and supplement Table 3 
beta_adj <- (t(beta_hat[ ,-1]) - alpha_s * beta) / (sigma_s / tau)

all_alpha <- c(0.99, 0.98, 0.95, 0.90, 0.8)
result <- matrix(0, 5, 2)
for(i in 1:5){
  alpha <- all_alpha[i]
  
  normal_quantile <- qnorm(0.5 + alpha / 2)
  covered <- beta_adj < normal_quantile & beta_adj > -normal_quantile
  
  result[i, ] <- c(mean(colMeans(covered)),
                   sd(colMeans(covered)) / sqrt(N)
  )
}

result * 100

# 3. p-value probability of a two-sided t-test
# Supplement Table 2 
j <- sample(which(beta == 0), 1) # sample a random null coordinate
p_val <- 2 * pnorm(abs(beta_adj[j, ]), lower.tail = F)
pval_probs <- c(0.1, 0.05, 0.01, 0.005)

result <- matrix(0, 4, 2)
for(i in 1:4){
  prob <- pval_probs[i]
  result[i, ] <- c(mean(p_val < prob),
                   sd(p_val < prob) / sqrt(N)
  )
}
result * 100



