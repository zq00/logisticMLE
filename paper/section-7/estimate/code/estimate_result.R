# Section 7.2. Estimating model parameters
# Results in Table 10 (II), Table 11 (II) and supplement Table 4 - 6
library(tidyverse)
library(glmhd)

# -- global theme for plots 
my_theme <- theme(axis.text.x = element_text(color = "grey10", size = 12),
                  axis.text.y = element_text(color = "grey10", size = 12),
                  axis.title.x = element_text(color = "black", size = 15),
                  axis.title.y = element_text(color = "black", size = 15),
                  plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm")) + theme_bw()
# --
n <- 4000
p <- 800
B <- 15
type = "ar5"

fileloc <- "/estimate/"  
filenames <- list.files(paste0(fileloc, "output/", type , "/"))

index <- sort(unique(as.numeric(unlist(sapply(filenames, function(t) strsplit(strsplit(t, "_")[[1]][2], ".txt"))))))
mle <- matrix(0, length(index) * B, p+1)
alphahat <- numeric(length(index) * B)
sigmahat <- matrix(0, length(index) * B, p)
# lrt <- numeric(length(index) * B)
for(i in 1:length(index)){
  new_name <- paste0(fileloc, "output/", type, "/mle_", index[i], ".txt")
  mle[((i-1)*B+1):(i*B), ] <- as.matrix(read.table(new_name))
  
  new_name <- paste0(fileloc, "output/", type, "/sigmahat_", index[i], ".txt")
  sigmahat[((i-1)*B+1):(i*B), ] <- as.matrix(read.table(new_name))
  
  new_name <- paste0(fileloc, "output/", type, "/alphahat_", index[i], ".txt")
  alphahat[((i-1)*B+1):(i*B)] <- scan(new_name)
}

beta <- scan(paste0(fileloc, "param/beta_", type,".txt"))
beta0 <- 1

B <- 10000 
mle <- mle[1:B, ]
alphahat <- alphahat[1:B]
sigmahat <- sigmahat[1:B, ]

# 1. coverage proportion of a non-null
# Table 10 Column (II), supplement Table 4
# Supplementary Table 5

j <- sample(which(beta != 0), 1) # randomly pick a non-null

# adjusted beta_hat
beta_adj <- (mle[,j+1] - alphahat * beta[j]) / sigmahat[,j]

# Coverage proportion
coverage_proportion <- tibble(
  prop = c(0.99, 0.98, 0.95, 0.90, 0.80),
  normal_quantile = qnorm(0.5 + prop/2)
) %>% 
  mutate(
    cov_prop = normal_quantile %>% map_dbl(~ mean(beta_adj < . & beta_adj > -.)),
    sd = normal_quantile %>% map_dbl(~ sd(beta_adj < . & beta_adj > -.) / sqrt(length(beta_adj)))
  ) 

coverage_proportion*100

# 2. p-value of a two-sided t-test
# Supplement Table 5
j <- sample(which(beta == 0), 1) # sample a null

beta_adj <- mle[,j+1] / sigmahat[,j]
p_val <- 2 * pnorm(-abs(beta_adj))
coverage_proportion <- tibble(
  prop = c(0.1, 0.05, 0.01, 0.005),
) %>% 
  mutate(
    cov_prop = prop %>% map_dbl(~ mean(p_val < . )),
    sd = prop %>% map_dbl(~ sd(p_val < . ) / sqrt(length(p_val)))
  ) 

coverage_proportion * 100

# 3.Coverage proportion of all the variables (bulk)
# Table 11 (Column II), supplement Table 6
beta_adj <- (t(mle[ ,-1]) - outer(beta, alphahat, "*" )) / t(sigmahat)

coverage_proportion <- tibble(
  prop = c(0.99, 0.98, 0.95, 0.90, 0.8),
  normal_quantile = qnorm(0.5 + prop/2)
) %>% 
  mutate(
    cov_prop = normal_quantile %>% map_dbl(~ mean(beta_adj < . & beta_adj > -.)),
    sd = normal_quantile %>% map_dbl(~ sd(colMeans(beta_adj < . & beta_adj > -.)) / sqrt(nrow(beta_adj)))
  ) 

coverage_proportion * 100

