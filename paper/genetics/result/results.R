#This script produces plot in Section 6 - genetic dataset
# run read_data.R first to load results
library(ggplot2) 
#-- parameters
kappa=1454/5000
#-- 

#-- global theme for plots
theme(axis.text.x = element_text(color = "grey10", size = 12),
      axis.text.y = element_text(color = "grey10", size = 12),
      axis.title.x = element_text(color = "black", size = 15),
      axis.title.y = element_text(color = "black", size = 15),
      plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm"))
#--
  
#-- obtain estimated parameters
setwd("/genetic/result")

B <- 5000
beta_est <- beta_est[1:B,]
cov_est <- cov_est[1:B,]
dev <- dev[1:B]
gamma_hat <- gamma_hat[1:B]
param <- param[1:B,]
pc <- pc[1:B,]
std_est <- std_est[1:B,]

hat_alpha = param[,1];hat_sigma = param[,2];hat_lambda=param[,3]
# --

j_null = 557

# -- Distribution of one null coordinate
#obtain hat(\beta)
hat_beta = beta_est[,j_null]
#obtain hat(\tau)
hat_tau = cov_est[,j_null]/sqrt(1-kappa)
#obtain hat(\sigma)
#compute scaled tilde(\beta)
tilde_beta = hat_beta * hat_tau / hat_sigma

#p-values from a two-sided t-test
p.t <- data.frame(2*pnorm(hat_tau*abs(hat_beta)/hat_sigma,0,1,lower.tail=FALSE))
colnames(p.t) <- c('x')
#figure 3 (b) histogram of p-value from t-test
g_pval_t <- ggplot(p.t) +
  geom_histogram(aes(x = x),
                 breaks = seq(0,1,by=0.05),fill='skyblue2', color='black',alpha=0.8)+
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey15", size = 12),
        axis.text.y = element_text(color = "grey15", size = 12),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm"))+
  ylab('Density')+xlab('P-values')

filename <- "/genetic/figure/pval_t.png"
ggsave(filename, plot = g_pval_t, device = "png",
       scale = 1, width = 4, height = 3, units = "in", 
       dpi = 300)

#figure 3 (c) empirical cdf for this estimated beta_hat
ecdf_beta = data.frame(cbind(sort(p.t$x,decreasing=FALSE),(seq(1,length(p.t$x),1)/(length(p.t$x)+1))))
colnames(ecdf_beta) = c('x','y')

theme_set(theme_bw())
g_ecdf_t <- ggplot(ecdf_beta, aes(x=x, y=y))+geom_point(size=0.8) + geom_abline(intercept=0,slope=1, linetype="dashed", color = "red")+
  theme(axis.text.x = element_text(color = "grey20", size = 12),
        axis.text.y = element_text(color = "grey20", size = 12),
        axis.title.x = element_text(color = "grey20", size = 15),
        axis.title.y = element_text(color = "grey20", size = 15),
        plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm"))+  ylab('Empirical cdf')+xlab('Sorted p-value')
g_ecdf_t

filename <- "/genetic/figure/ecdf_t.png"
ggsave(filename, plot = g_ecdf_t, device = "png",
       scale = 1, width = 4, height = 3, units = "in", 
       dpi = 300)

#table 7 column 1
mean(p.t<=0.1)*100; sqrt(mean(p.t<=0.1)*(1-mean(p.t<=0.1))/B)*100
mean(p.t<=0.05)*100;sqrt(mean(p.t<=0.05)*(1-mean(p.t<=0.05))/B)*100
mean(p.t<=0.01)*100;sqrt(mean(p.t<=0.01)*(1-mean(p.t<=0.01))/B)*100
mean(p.t<=0.001)*100;sqrt(mean(p.t<=0.001)*(1-mean(p.t<=0.001))/B)*100

# 2-sided t-test from R
#figure 3 (a): p-value from t-test from R
p.tr <- data.frame(2*pnorm(abs(hat_beta)/std_est[,j_null],0,1,lower.tail=FALSE))
colnames(p.tr) <- c('x')

g_pval_tr <- ggplot(p.tr) +
  geom_histogram(aes(x = x),
                 breaks = seq(0,1,by=0.05),fill='skyblue2', color='black',alpha=0.8)+
  theme(axis.text.x = element_text(color = "grey15", size = 12),
        axis.text.y = element_text(color = "grey15", size = 12),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm"))+
  ylab('Density')+xlab('P-values')

filename <- "/genetic/figure/pval_tr.png"
ggsave(filename, plot = g_pval_tr, device = "png",
       scale = 1, width = 4, height = 3, units = "in", 
       dpi = 300)

#-- Distribution of the LLR for the same null coordinate
#obtain llr
#obtain hat(\lambda)
#obtain p-value according to kappa * hat(\sigma)^2 / hat(\lambda) * chi2_1
p.val = pchisq(dev * hat_lambda / kappa / hat_sigma^2, df = 1, lower.tail=FALSE)

# table 7 column 3
mean(p.val<=0.1)*100; sqrt(mean(p.val<=0.1)*(1-mean(p.val<=0.1))/B)*100
mean(p.val<=0.05)*100;sqrt(mean(p.val<=0.05)*(1-mean(p.val<=0.05))/B)*100
mean(p.val<=0.01)*100;sqrt(mean(p.val<=0.01)*(1-mean(p.val<=0.01))/B)*100
mean(p.val<=0.001)*100;sqrt(mean(p.val<=0.001)*(1-mean(p.val<=0.001))/B)*100

#figure 4 (c):empirical cdf of adjusted p-values
ecdf_p = data.frame(cbind(sort(p.val,decreasing=FALSE),(seq(1,length(p.val),1)/(length(p.val)+1))))
colnames(ecdf_p) = c('x','y')

theme_set(theme_bw())
g_ecdf_llr <- ggplot(ecdf_p, aes(x=x, y=y))+geom_point(size=0.8) + 
  geom_abline(intercept=0,slope=1, linetype="dashed", color = "red")+
  theme(axis.text.x = element_text(color = "grey20", size = 12),
        axis.text.y = element_text(color = "grey20", size = 12),
        axis.title.x = element_text(color = "grey20", size = 15),
        axis.title.y = element_text(color = "grey20", size = 15),
        plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm"))+  ylab('Empirical cdf')+xlab('P-values')
plot(g_ecdf_llr)

filename <- "/genetic/figure/ecdf_llr.png"
ggsave(filename, plot = g_ecdf_llr, device = "png",
       scale = 1, width = 4, height = 3, units = "in", 
       dpi = 300)

#figure 6 (b) : histogram of the p-values for the llr
p.val <- data.frame(p.val)
colnames(p.val) <- c('x')

g_pval_llr <- ggplot(p.val) +
  geom_histogram(aes(x = x),
                 breaks = seq(0,1,by=0.05),fill='skyblue2', color='black',alpha=0.8)+
  theme(axis.text.x = element_text(color = "grey15", size = 12),
        axis.text.y = element_text(color = "grey15", size = 12),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm"))+
  ylab('Density')+xlab('P-values')
  
filename <- "/genetic/figure/g_pval_llr.png"
ggsave(filename, plot = g_pval_llr, device = "png",
       scale = 1, width = 4, height = 3, units = "in", 
       dpi = 300)

#figure 6 (a) :  histogram of the unadjusted p-values
p.unj = pchisq(dev, df = 1, lower.tail=FALSE)
p.unj <- data.frame(p.unj)
colnames(p.unj) <- c('x')

g_pval_llr_unadj <- ggplot(p.unj) +
  geom_histogram(aes(x = x),
                 breaks = seq(0,1,by=0.05),fill='skyblue2', color='black',alpha=0.8)+
  theme(axis.text.x = element_text(color = "grey15", size = 12),
        axis.text.y = element_text(color = "grey15", size = 12),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm"))+
  ylab('Density')+xlab('P-values')

filename <- "/genetic/figure/pval_llr_unadj.png"

ggsave(filename, plot = g_pval_llr_unadj, device = "png",
       scale = 1, width = 4, height = 3, units = "in", 
       dpi = 300)

#-- Coverage
# Table 8: coverage proportion of one non-null variable
beta <- scan("/genetics/param/beta.txt")
non_null <- scan("/genetics/param/non_null.txt")
#obtain hat(\beta)
hat_beta = as.matrix(beta_est)
#obtain hat(\tau)
hat_tau = as.matrix(cov_est)/sqrt(1-kappa)

tilde_beta <- matrix(0, B, 1454)
for(i in 1:B){
  tilde_beta[i,] <- hat_tau[i,] * (hat_beta[i,]-hat_alpha[i]*beta)/hat_sigma[i] 
}

# randomly pick a non-null
j_non_null <- sample(non_null, 1)

p_0.99 <- mean(abs(tilde_beta[,j_non_null])<2.58); 
p_0.99*100; sqrt(p_0.99*(1-p_0.99)/B)*100 #99%
p_0.98 <- mean(abs(tilde_beta[,j_non_null])<2.32); 
p_0.98*100; sqrt(p_0.98*(1-p_0.98)/B)*100 #98%
p_0.95 <- mean(abs(tilde_beta[,j_non_null])<1.96); 
p_0.95*100; sqrt(p_0.95*(1-p_0.95)/B)*100 #95%
p_0.90 <- mean(abs(tilde_beta[,j_non_null])<1.64); 
p_0.90*100; sqrt(p_0.90*(1-p_0.90)/B)*100 #90%
p_0.80 <- mean(abs(tilde_beta[,j_non_null])<1.28); 
p_0.80*100;sqrt(p_0.80*(1-p_0.80)/B)*100 #80%

# Figure 5
betahat_nonnull <-  data.frame(
  cbind(
    sort(tilde_beta[,j_non_null],decreasing=FALSE),
    qnorm((seq(1,B,1)/(B+1))),0,1)
  )
colnames(betahat_nonnull) = c('x','y')

g_qqplot <- ggplot(betahat_nonnull, aes(x=x, y=y))+
  geom_point(size=0.8) + 
  geom_abline(intercept=0,slope=1, linetype="dashed", color = "red")+
  theme(axis.text.x = element_text(color = "grey20", size = 12),
        axis.text.y = element_text(color = "grey20", size = 12),
        axis.title.x = element_text(color = "grey20", size = 15),
        axis.title.y = element_text(color = "grey20", size = 15),
        plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm"))+  
  ylab('Empirical quantiles')+
  xlab('Normal quantiles')
plot(g_qqplot)

filename <- "/genetics/figure/subgaussian-nonnull.png"
ggsave(filename, plot = g_qqplot, device = "png",
       scale = 1, width = 4, height = 3, units = "in", 
       dpi = 300)

#coverage proportion overall
#coverage proportion
mean(rowMeans(abs(tilde_beta[,non_null])<2.32))*100; sd(rowMeans(abs(tilde_beta[,non_null])<2.32))/sqrt(B)*100 #98%
mean(rowMeans(abs(tilde_beta[,non_null])<1.96))*100; sd(rowMeans(abs(tilde_beta[,non_null])<1.96))/sqrt(B)*100 #95%
mean(rowMeans(abs(tilde_beta[,non_null])<1.64))*100; sd(rowMeans(abs(tilde_beta[,non_null])<1.64))/sqrt(B)*100 #90%
mean(rowMeans(abs(tilde_beta[,non_null])<1.28))*100; sd(rowMeans(abs(tilde_beta[,non_null])<1.28))/sqrt(B)*100 #80%

null <- (1:1454)[-non_null]
mean(rowMeans(abs(tilde_beta[,null])<2.32))*100; sd(rowMeans(abs(tilde_beta[,null])<2.32))/sqrt(B)*100 #98%
mean(rowMeans(abs(tilde_beta[,null])<1.96))*100; sd(rowMeans(abs(tilde_beta[,null])<1.96))/sqrt(B)*100 #95%
mean(rowMeans(abs(tilde_beta[,null])<1.64))*100; sd(rowMeans(abs(tilde_beta[,null])<1.64))/sqrt(B)*100 #90%
mean(rowMeans(abs(tilde_beta[,null])<1.28))*100; sd(rowMeans(abs(tilde_beta[,null])<1.28))/sqrt(B)*100 #80%

mean(rowMeans(abs(tilde_beta)<2.32))*100; sd(rowMeans(abs(tilde_beta)<2.32))/sqrt(B)*100 #98%
mean(rowMeans(abs(tilde_beta)<1.96))*100; sd(rowMeans(abs(tilde_beta)<1.96))/sqrt(B)*100 #95%
mean(rowMeans(abs(tilde_beta)<1.64))*100; sd(rowMeans(abs(tilde_beta)<1.64))/sqrt(B)*100 #90%
mean(rowMeans(abs(tilde_beta)<1.28))*100; sd(rowMeans(abs(tilde_beta)<1.28))/sqrt(B)*100 #80%

# -- Below are some other plots 
# histogram of coverage proportion 
pcov <- data.frame(pc)
colnames(pcov) <- c('x','y')

ggplot(pcov) +
  geom_histogram(aes(x = y),
                 breaks = seq(0.8,1,by=0.005),fill='skyblue2', color='black',alpha=0.8)+
  geom_histogram(aes(x = x),
                 breaks = seq(0.75,1,by=0.005),fill='lightpink2', color='black',alpha=0.8)+
  geom_segment(aes(x = 0.95, y = 0, xend = 0.95, yend = 1300), color="red", linetype="dashed")+
  theme(axis.text.x = element_text(color = "grey15", size = 12),
        axis.text.y = element_text(color = "grey15", size = 12),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm"))+
  scale_y_continuous(breaks=seq(0, 1250, 250))+
  ylab('Counts')+xlab('Coverage proportion')



