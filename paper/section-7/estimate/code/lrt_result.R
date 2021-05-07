# Section 7.2 Estimating model parameters - LRT
# Table 12
library(glmhd)
library(tidyverse)

# read the deviance statistics and the estimated parameters
fileloc <- paste0("/estimate/output/ar5/")
filenames <- list.files(fileloc) 
filenames <- filenames[grepl("lrt", filenames)]
indices <- sapply(filenames, function(t) strsplit(t, "_")[[1]][2])

dev <- paste0("lrt_", indices) %>% 
  map(~ read.table(paste0(fileloc, .))) %>% 
  bind_rows()
params <- paste0("param_", indices) %>% 
  map(~ read.table(paste0(fileloc, .))) %>% 
  bind_rows()

# use a total of 10,000 simulations
B <- 10000
dev <- dev[1:B, ]
params <- params[1:B, ]

# 1. theoretical adjustment
# Table 12 (Column I)
# param <- find_param(kappa = 0.2, gamma = sqrt(5), intercept = T, beta0 = 1)
adj <- param[2] / param[3]^2  # theoretical adjustment
# compute the adjusted p-values
padj <- pchisq(dev * adj, 1, lower.tail = F)

pvalProb <- tibble(
  prop = c(0.1, 0.05, 0.01, 0.005),
) %>% 
  mutate(
    probs = prop %>% map_dbl(~ mean(padj<.)),
    sd = probs %>% map_dbl(~ sd(padj<.) / sqrt(length(padj)))
  ) 

pvalProb*100

# 2. Estimated parameters
# Table 12 (Column II)
adj <- params[,2] / params[,3]^2  # estimated adjustment factors

# compute the adjusted p-values
padj <- pchisq(dev * adj, 1, lower.tail = F)

pvalProb <- tibble(
  prop = c(0.1, 0.05, 0.01, 0.005),
) %>% 
  mutate(
    probs = prop %>% map_dbl(~ mean(padj<.)),
    sd = probs %>% map_dbl(~ sd(padj<.) / sqrt(length(padj)))
  ) 

pvalProb*100

# 3. Classical (no adjustment)
# Table 12 (Column III)

pval <- pchisq(dev, 1, lower.tail = F)

pvalProb <- tibble(
  prop = c(0.1, 0.05, 0.01, 0.005),
) %>% 
  mutate(
    probs = prop %>% map_dbl(~ mean(pval<.)),
    sd = probs %>% map_dbl(~ sd(pval<.) / sqrt(length(pval)))
  ) 

pvalProb*100

