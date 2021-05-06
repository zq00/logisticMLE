# Section 3.2.2  Condition on the regression coef
# Simulation results (Figure 4)
library(glmhd)
library(tidyverse)

# -- global theme for plots 
my_theme <- theme(axis.text.x = element_text(color = "grey10", size = 12),
                  axis.text.y = element_text(color = "grey10", size = 12),
                  axis.title.x = element_text(color = "black", size = 15),
                  axis.title.y = element_text(color = "black", size = 15),
                  plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm")) + theme_bw()
# --

kappa <- 0.2
gamma2 <- 4
n <- c(1000,2000, 4000,8000)

fileloc <- "/output/"
# for each parameter sequence 
# theoretical parameter 
params <- find_param(kappa = kappa, gamma = sqrt(gamma2), beta0 = 0, intercept = F)

for(nn in n){
  alpha <- NULL
  folder_name <- paste(kappa, gamma2, nn, sep = "-")
  all_files <- list.files(paste0(fileloc, folder_name, "/"))
  for(i in 1:length(all_files)){
    new_file <- scan(paste0(fileloc,folder_name,"/", all_files[i]), quiet = T)
    alpha <- c(alpha, new_file)
  }
  # compute the standard deviation
  cat("_____")
  cat("kappa= ", kappa, ", gamma = ", sqrt(gamma2), " n = ", nn, "\n")
  cat("Number of simulations is", length(alpha), "\n")
  cat("Observed bias = ",sqrt(nn) * abs(mean(alpha) - params[1]), "\n" )
  cat("observed std = ", sqrt(nn) * sqrt(mean((alpha - params[1])^2)), "\n")
  cat("observed range of alpha(n) is ", round(sqrt(nn) * range(alpha - params[1]), 3), "\n")
}


# Plot the standard deviations
data <- tibble(
  group = factor(rep(c(1,2,3), each = 4)),
  kappa = rep(c(0.1, 0.2), time = c(8, 4)),
  gamma = rep(c(2, 5, 2), each = 4),
  n = rep(c(1, 2, 4, 8), time = 3),
  sd = c(2.80, 2.77, 2.74, 2.77, 
         5.62, 5.30, 5.21, 5.11,
         4.77, 4.63, 4.62, 4.55)
)

filename <- "/result/result.png"
g <- ggplot(data = data) + 
  geom_point(aes(x = n, y = sd, color = group)) + 
  geom_line(aes(x = n, y = sd, color = group)) +
  ylab(expression(paste(hat(sd),"(", sqrt(n)(alpha[n] - alpha[s]), ")"))) +
  scale_x_continuous("N (/1000)", breaks = c(1,2,4,8), labels = c("1","2","4","8")) +
  scale_color_discrete(name = "Parameters", labels = c(expression(paste(kappa, "=0.1", ", ",gamma,"=2")),
                                                       expression(paste(kappa, "=0.1", ", ",gamma,"=5")),
                                                       expression(paste(kappa, "=0.2", ", ",gamma,"=2")))) + 
  my_theme
ggsave(filename, plot = g, device = "png",
       scale = 1, width = 5, height = 3, units = "in", 
       dpi = 300)
