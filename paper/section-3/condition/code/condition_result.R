# Section 3.2.2 Condition on the regression coef
# code to generate Figures from simulated MLE
library(tidyverse)
# -- global theme for plots 
my_theme <- theme(axis.text.x = element_text(color = "grey10", size = 12),
                  axis.text.y = element_text(color = "grey10", size = 12),
                  axis.title.x = element_text(color = "black", size = 15),
                  axis.title.y = element_text(color = "black", size = 15),
                  plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm")) + theme_bw()
# --

params <- data.frame(
  ratio = rep(c(0.15, 0.3, 0.5), each = 4),
  n = rep(c(1000, 2000, 4000, 8000), 3),
  nonnull = c(80, 106, 133, 93,
              195, 30, 592, 825,
              3, 162, 793, 1576)
)

betaloc <- "/Users/zq/Documents/GitHub/logisticMLE/paper/section-3/condition/param"
fileloc <- "/Users/zq/Documents/Simulation_Data/glm/condition4/output/"
gamma2 <- 4
kappa <- 0.2
# for each parameter sequence
for(k in 1:12){
  n <- params[k, 2]
  r <- params[k, 1]
  p <- n * kappa
  nonnull <- params[k, 3]
  
  folder_name <- paste(n, r, sep = "-")
  all_files <- list.files(paste0(fileloc, folder_name, "/"))
  mle <- NULL
  if(length(all_files) == 0){
    next;
  }else{
    # load results
    for(i in 1:length(all_files)){
      new_file <- scan(paste0(fileloc,folder_name,"/",all_files[i]), quiet = T)
      mle <- c(mle, new_file)
    }
    # load beta
    filename <- paste0(paste(n, r, sep = "_"), ".txt")
    beta <- scan(paste0(betaloc, filename))
    
    cat("______________")
    cat("ratio = ", r, ",  n = ", n, "\n")
    cat("n simulation = ", length(mle), "\n")
    cat("observed bias is = ", mean(mle) / beta[nonnull], "\n")
    cat("Std is", round(sd(mle) / sqrt(kappa), 3), "\n")
    cat("fourth moment is", (round(mean((mle-mean(mle))^4)))/sqrt(kappa), "\n")
    # normal quantile plot of the MLE of a nonnull coefficient
    g <-  ggplot(data = tibble(mle = mle), aes(sample = mle))+
      stat_qq(size = 0.8) + 
      stat_qq_line(linetype="dashed", color = "red") + 
      xlab('Normal quantiles')+ylab('Empirical quantiles')+
      my_theme
    
    filename <- paste0("/Users/zq/Documents/Simulation_Data/glm/condition4/result/",folder_name, ".png") # file location
    ggsave(filename, plot = g, device = "png",
           scale = 1, width = 4, height = 3, units = "in", 
           dpi = 300)
  }
}

# Figure 3 (a)
# plot standard deviations versus n for different ratio 
data <- tibble(
  group = factor(rep(c(1,2,3), each = 4)),
  ratio = rep(c(0.15, 0.3, 0.5), each = 4),
  n = rep(c(1, 2, 4, 8), time = 3),
  sd = c(6.01, 5.96, 5.90, 5.85,
         6.59, 6.49, 6.49, 6.42, 
         7.83, 7.76, 7.71, 7.75)
)

# params <- find_param(kappa = 0.2, gamma = 2, intercept = F) # compute theoretical parameter
sigma_theory <- 1.956293 / sqrt(0.2) / 0.77460
filename <- "/Users/zq/Documents/Simulation_Data/glm/condition4/result/result_incomplete.png"

g <- ggplot(data = data) + 
  geom_point(aes(x = n, y = sd, color = group)) + 
  geom_line(aes(x = n, y = sd, color = group)) +
  geom_abline(slope = 0, intercept = sigma_theory, color = "black", linetype = "dashed") + 
  ylim(c(5.5, 8)) +
  ylab(expression(paste("SD(",hat(beta[j]),")"))) +
  scale_x_continuous("N (/1000)", breaks = c(1,2,4,8), labels = c("1","2","4","8")) +
  scale_color_discrete(name = expression(paste(tau[j], beta[j], "/", gamma)), labels = c(0.15, 0.3, 0.5)) + 
  my_theme
g
ggsave(filename, plot = g, device = "png",
       scale = 1, width = 5, height = 3, units = "in", 
       dpi = 300)




