library(tidyverse)

#---Global theme for figures 
my_theme <- theme(axis.text.x = element_text(color = "grey10", size = 12),
      axis.text.y = element_text(color = "grey10", size = 12),
      axis.title.x = element_text(color = "black", size = 15),
      axis.title.y = element_text(color = "black", size = 15),
      plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm")) + theme_bw()
#---

type = "ar_0.5"
fileloc <- paste0("/section2.3/output/", type , "/")
filenames <- list.files(fileloc)

# Estimated beta_hats
beta_hat <- filenames %>% 
  map(~ read.table(paste0(fileloc, .))) %>% 
  bind_rows()

# Parameters
n <- 4000
p <- 800
gamma2 <- 5

# True coefficients
Sigma <- as.matrix(read.table(
  paste0("/section2.3/param/", type, "/Sigma.txt")))
non_null <- scan(paste0("/section2.3/param/", type, "/non_null.txt"))

beta <- numeric(p)
beta_solve <- uniroot(function(x){ beta[non_null] <- x; return(t(beta)%*%(Sigma%*%beta)/n - gamma2);}, c(0,300), check.conv=TRUE,extendInt = "upX")$root
beta[non_null] <- beta_solve

j_nonnull <- scan(paste0("/section2.3/param/", type, "/j_nonnull.txt"))
beta <- beta[non_null[1]]

# tau_j
tau <- sqrt(1/diag(solve(Sigma)))[j_nonnull]

# empirical distribution
alpha_s = 1.4994
sigma_s = 4.7436

# adjusted beta_hat
beta_adj <- (beta_hat[,2] - alpha_s * beta) * tau[2] / sigma_s

# normal quantile plot 
data <- tibble(
  index = 1:100000,
  beta_hat = sort(beta_adj[1:100000])
) %>% 
  mutate(
    quantile = index / 100000,
    normal_quantile = qnorm(quantile, 0, 1)
  ) 

g <- data %>% 
  ggplot(aes(x=normal_quantile, y=beta_hat))+
  geom_point(size=0.8) + 
  geom_abline(intercept=0,slope=1, linetype="dashed", color = "red")+
  xlab('Normal quantiles')+ylab('Empirical quantiles')+
  my_theme
plot(g)

filename <- "/section2.3/Figure/sec2_3.png"
ggsave(filename, plot = g, device = "png",
       scale = 1, width = 4, height = 3, units = "in", 
       dpi = 300)

# Coverage proportion
coverage_proportion <- tibble(
  prop = c(0.99, 0.98, 0.95, 0.90, 0.8),
  normal_quantile = qnorm(0.5 + prop/2)
) %>% 
  mutate(
    cov_prop = normal_quantile %>% map_dbl(~ mean(beta_adj < . & beta_adj > -.)),
    sd = normal_quantile %>% map_dbl(~ sd(beta_adj < . & beta_adj > -.) / sqrt(100000))
) 

coverage_proportion*100



  




