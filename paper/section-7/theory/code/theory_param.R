# Section 7.1.2. Effect of the intercept
# Compute the theoretical parameters (Table 13)
library(glmhd)

# Column I
find_param(kappa = 0.2, gamma = sqrt(5), beta0 = 0, intercept = T)
find_param(kappa = 0.2, gamma = sqrt(5), beta0 = 0.5, intercept = T)
find_param(kappa = 0.2, gamma = sqrt(5), beta0 = 1, intercept = T)
find_param(kappa = 0.2, gamma = sqrt(5), beta0 = 2, intercept = T)
find_param(kappa = 0.2, gamma = sqrt(5), beta0 = 2.5, intercept = T)

# Column II
find_param(kappa = 0.2, gamma = sqrt(5), intercept = F)
find_param(kappa = 0.2, gamma = sqrt(5 + 0.5^2), intercept = F)
find_param(kappa = 0.2, gamma = sqrt(5 + 1), intercept = F)
find_param(kappa = 0.2, gamma = sqrt(5 + 4), intercept = F)
find_param(kappa = 0.2, gamma = sqrt(5 + 2.5^2), intercept = F)