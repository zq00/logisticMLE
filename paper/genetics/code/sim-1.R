#genetic example
library(SNPknock)
load('/genetics/param/hmm.RData')

# -- Job submission arguments
args = commandArgs(trailingOnly=TRUE)
ind <- args[1]

sim_ind = 8;
# -- 

fileloc = paste0("/genetic/output_",sim_ind,"/data/data_",ind,"/")
print(fileloc)
system(paste0('mkdir -p ',fileloc));

phase0 = read.table("/section-5/code/phase0.txt")

# -- Parameters
n <- 5000;
kappa = 1454/n; p =1454;

# -- Parameters ProbeFrontier
B = 10; #num of repetition
# --

beta = scan('/genetic/param/beta.txt')
non_null = scan('/genetic/param/non_null.txt')

# samples one observation
temp = floor(sim_ind*1000+as.numeric(ind)); print(temp) #debug
set.seed(temp)
XX = SNPknock.models.sampleHMM(hmm$pInit, hmm$Q, hmm$pEmit,n=n)

X = apply(XX, 2, function(x) {return((x-mean(x))/sd(x))})/sqrt(n)
print(X[1:2,1:9]) #debug
eta <- X %*% beta; mu <- exp(eta) / (1+exp(eta))
print(mean(eta^2)); #check what is gamma^2
y <- rbinom(n, 1, mu)

#write data             
write.table(X,paste0(fileloc,"X.txt"),col.names = FALSE, row.names=FALSE);
write.table(y,paste0(fileloc,"Y.txt"),col.names = FALSE, row.names=FALSE);

# ProbeFrontier estimate of gamma2
kappas = seq(kappa, 0.5, by = 0.01); n_kappa = length(kappas);
pi_hats = numeric(n_kappa);
is_sep = numeric(B)

# -- Begin ProbeFrontier
for(i in 1:n_kappa){
	stop=FALSE;
	for(b in 1:B){
		k = kappas[i]; nn = floor(p / k); #number of samples in subsample

		sind = sample(1:n, nn, replace = FALSE); #subsample without replacement
		Xs = X[sind,];Ys = y[sind];
		write.table(Xs,paste0(fileloc,'Xs.txt'),col.names = FALSE, row.names=FALSE);
		write.table(Ys,paste0(fileloc,'Ys.txt'),col.names = FALSE, row.names=FALSE);

		#check if mle exists
		command = paste0('ml load matlab; matlab -nodisplay -r "addpath /scratch/users/qzhao1/logistic/ProbeFrontier/script; fileloc = \'', fileloc, '\';is_sep; exit;"')
		system(command)
		#obtain result
		is_sep[b] = scan(paste0(fileloc,'sep.txt'))
		print(c(b,is_sep[1:b]));
	}
	pi_hats[i] = mean(is_sep);
	print(c(k,pi_hats[1:i]));
	if(pi_hats[i]>=0.5){
		stop = TRUE; #put a flag that have found kappa                  

		kappa_n = seq(kappas[i]-0.01, kappas[i], by=0.001) #check more values
		pi_n = numeric(length(kappa_n))

		#nov.18
		#better control
		pi_n[1] = pi_hats[i-1]; pi_n[length(kappa_n)] = pi_hats[i]
		kappa_hat = kappas[i];

		for(j in 2:(length(kappa_n)-1)){
			k = kappa_n[j]; nn =floor(p/k);
			for(b in 1:B){
				sind = sample(1:n, nn, replace = FALSE);
				Xs = X[sind,];Ys = y[sind];
				write.table(Xs,paste0(fileloc,'Xs.txt'),col.names = FALSE, row.names=FALSE);
				write.table(Ys,paste0(fileloc,'Ys.txt'),col.names = FALSE, row.names=FALSE);
				#check if mle exists
				command = paste0('ml load matlab; matlab -nodisplay -r "addpath /section-5/code/; fileloc = \'', fileloc, '\';is_sep; exit;"')
				system(command)
				#obtain result
				is_sep[b] = scan(paste0(fileloc,'sep.txt'))
				print(c(k,is_sep[1:b]));
			}
			pi_n[j] = mean(is_sep)
			print(c(pi_n[1:j]));
		if(pi_n[j]>=0.5){
				kappa_hat = kappa_n[j-1] + (0.5 - pi_n[j-1])/(pi_n[j] - pi_n[j-1]) * (kappa_n[j] - kappa_n[j-1]);
				#kappa_hat = kappa_n[j];
				print(paste0("have found kappa_hat ",kappa_hat));break;
			}
		}
	}
	if(stop==TRUE){break;}
}
print(kappa_hat)
gamma_hat = min(phase0[which(phase0[,2]<=kappa_hat),1]);
print(paste0("estimated gamma is", gamma_hat));

write.table(gamma_hat,paste0(fileloc,'gamma.txt'),col.names=FALSE,row.names=FALSE);

command = paste0('ml load matlab; matlab -nodisplay -r "addpath /section-5/code/; fileloc = \'', fileloc, '\';kappa=',kappa,'; FindParam; exit;"')
system(command);

param_new = as.matrix(read.table(paste0(fileloc,"param.txt")));

print(param);
# -- End ProbeFrontier

# -- Run Logistic regression on (X,y)
cov_store <- 1/sqrt(diag(solve(t(X) %*% X)))

# fit a logistic regression
fit <- glm(y~X+0, family=binomial)
beta_est <- fit$coefficients
std_est <- summary(fit)$coef[,2]

fit0 <- glm(y~X[,-j_null]+0, family=binomial)
dev_store<- -anova(fit, fit0)$Deviance[2]

#compute coverage proportion
ci_r <- cbind( fit$coefficients - 1.96 * std_est,fit$coefficients +  1.96 * std_est)
pc_r <- mean((ci_r[,2] >= beta) & (ci_r[,1] <= beta))

print(c(pc_est,pc_r))

#ProbeFrontier estimaties
alpha_hat = param_new[1,1]
sigma_hat = param_new[1,2]
lambda_hat = param_new[1,3]

ci_pf = cbind(fit$coefficients - 1.96 * sigma_hat/(cov_store/(sqrt(1-kappa))), fit$coefficients + 1.96 * sigma_hat / (cov_store/sqrt(1-kappa))) / alpha_hat
pc_pf = mean((ci_pf[,2]>=beta) & (ci_pf[,1]<=beta))

print(paste0("simulation",bb,"coverage by ProbeFrontier is", pc_pf));




