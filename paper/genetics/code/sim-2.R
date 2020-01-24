# Fits MLE based on estimated parameters from the ProbeFrontier algorithm
# sim-1.R needs to be run before this

#specify file locations (please ignore)
sim_ind = 5;
fileloc = paste0('/genetic/output_',sim_ind,'/data/')

# -- parameters
n <- 5000
p = 1454
kappa = p/n

beta = scan('/genetic/param/beta.txt')
non_null = scan('/genetic/param/non_null.txt')

j_null = 557

# store estimated values 
# extracts all the successfully executed simulation indices 
B = scan(paste0('/genetics/exists_ex_',sim_ind,'.txt'))
nB = length(B)

beta_est = matrix(0,nB, p)
std_est = matrix(0,nB,p)
cov_est = matrix(0, nB, p)

dev_store = numeric(nB)
param = matrix(0,nB,3)
gamma_hat = numeric(nB)

pc = matrix(0,nB,2)

ci_r = NULL;ci_pf  =NULL #initialize to avoid memory issue

ii = 1

for(ind in B){
	#scan through all the indices
	newfile = paste0(fileloc, 'data_',ind,'/')
  
  	X = as.matrix(read.table(paste0(newfile,'X.txt')))
  	y = scan(paste0(newfile,'Y.txt'))

  	cov_est[ii,] <- 1/sqrt(diag(solve(t(X) %*% X)))
  
  	#fit a logistic regression
  	fit <- glm(y~X+0, family=binomial)
  	beta_est[ii,] <- fit$coefficients
  	std_est[ii,] <- summary(fit)$coef[,2]
	
  	fit0 <- glm(y~X[,-j_null]+0, family=binomial)
  	dev_store[ii] <- -anova(fit, fit0)$Deviance[2]
	ci_r <- cbind( fit$coefficients - 1.96 * std_est[ii,],fit$coefficients +  1.96 * std_est[ii,])
		
	#compute coverage proportion
	pc[ii,1] <- mean((ci_r[,2] >= beta) & (ci_r[,1] <= beta)) #1-r
		
	param_new = as.matrix(read.table(paste0(newfile,"param.txt")))	
		
	alpha_hat = param_new[1,1]
	sigma_hat = param_new[1,2]
	lambda_hat = param_new[1,3]
  
  	# CI from probefrontier estimates
	ci_pf = cbind(fit$coefficients - 1.96 * sigma_hat/(cov_est[ii,]/(sqrt(1-kappa))), fit$coefficients + 1.96 * sigma_hat / (cov_est[ii,]/sqrt(1-kappa))) / alpha_hat
	pc[ii,2] = mean((ci_pf[,2]>=beta) & (ci_pf[,1]<=beta)) #2-ProbeFrontier

 	print(paste0('ii=',ii,',alpha_hat=',alpha_hat)) #debug	    	
  	print(paste0('pc_pf=',pc[ii,2],', pc_r=',pc[ii,1])) #debug	
  	param[ii,] = param_new[1,]
  	gamma_hat[ii] = scan(paste0(newfile,'gamma.txt'))
  
  	ii = ii+1;
}

#store files 
fileloc = paste0('/genetic/output_',sim_ind,'/output/')

write.table(beta_est, paste0(fileloc,'beta_est_',sim_ind,'.txt'),col.names=FALSE, row.names=FALSE,append=FALSE)
write.table(cov_est, paste0(fileloc,'cov_est_',sim_ind,'.txt'),col.names=FALSE, row.names=FALSE,append=FALSE)
write.table(std_est, paste0(fileloc,'std_est_',sim_ind,'.txt'),col.names=FALSE, row.names=FALSE,append=FALSE)
write.table(dev_store, paste0(fileloc,'dev_',sim_ind,'.txt'),col.names=FALSE, row.names=FALSE,append=FALSE)
write.table(pc, paste0(fileloc,'pc_',sim_ind,'.txt'),col.names=FALSE, row.names=FALSE,append=FALSE)
write.table(param, paste0(fileloc,'param_',sim_ind,'.txt'),col.names=FALSE, row.names=FALSE,append=FALSE)
write.table(gamma_hat, paste0(fileloc,'gamma_hat_',sim_ind,'.txt'),col.names=FALSE, row.names=FALSE,append=FALSE)



