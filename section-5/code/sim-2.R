#find parameters: param are for all the parameters in 1-1000 etc.
#calculate beta_hat for all the data that we have estimated parameters

#-- Job submission arguments
args = commandArgs(trailingOnly=TRUE)
ii <- as.numeric(args[1])
B = as.numeric(args[2])
print(ii)
#-- 

fileloc = paste0('/section-5/output/',ii,'/') # store estimated beta_hat

# pick out the sucessfully completed jobs 
notexist = NULL
for(ind in 1:B){
  #scan through all the indices
  newfile = paste0(fileloc,'data_' ,ind,'/')
  param_new = try(as.matrix(read.table(paste0(newfile,"param.txt"))));
  if(class(param_new) == "try-error"){notexist=c(notexist,ind);}
}

exists= (1:B)[-notexist]
# --

# -- find estimated gamma_hat
param_est = NULL
gamma_hat = NULL

# -- calculate beta_hat for all the data that we have estimated parameters
indices = exists
B1 = length(indices)

kappa = 0.2
p=800

j_null = scan('/section-5/param/jnull.txt') 

beta_est_1 = matrix(0,B1 ,p)
cov_est_1 = matrix(0, B1,p)
dev_1 = numeric(B1)
std_est_1 = matrix(0, B1,p)
rho_hat = numeric(B1)

for(ind in 1:B1){
	b = indices[ind]
	print(c(ind,b))
	#read data
	newloc= paste0('/section-5/data/data_',b,'/')
	X = as.matrix(read.table(paste0(newloc, 'X.txt')))
	y = scan(paste0(newloc,'Y.txt'))
	#fit a logistic regression
 	fit <- glm(y~X+0, family=binomial)

 	beta_new <- fit$coefficients; beta_est_1[ind,] = beta_new
  	std_new <- summary(fit)$coef[,2]; std_est_1[ind,] = std_new

  	fit0 <- glm(y~X[,-j_null]+0, family=binomial)
  	dev_1[ind] <- -anova(fit, fit0)$Deviance[2]
    cov_new <- 1/sqrt(diag(solve(t(X) %*% X)))/(sqrt(1-kappa)); cov_est_1[ind,] = cov_new

    #estimate rho from the MLE
    rho_hat[ind] =sum(rowSums(X[,2:p]*X[,1:(p-1)]))/sum(rowSums(X[,2:p]^2))
    print(c(ind,rho_hat[ind])) 
    
    param_new = as.matrix(read.table(paste0(newloc,"param.txt")))
    param_est = rbind(param_est, param_new)
    
    gamma_new = as.matrix(read.table(paste0(newloc,"gamma.txt")))
    gamma_hat = c(gamma_hat, gamma_new)
}

fname = '/section-5/output/'
write.table(beta_est_1, paste0(fname,'beta_est_',ii,'.txt'),col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(std_est_1, paste0(fname,'std_est_',ii,'.txt'),col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(cov_est_1, paste0(fname,'cov_est_',ii,'.txt'),col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(dev_1, paste0(fname,'dev_',ii,'.txt'),col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(rho_hat, paste0(fname,'rho_',ii,'.txt'),col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(param_est, paste0(fname,'param_',ii,'.txt'),col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(gamma_hat, paste0(fname,'gamma_',ii,'.txt'),col.names=FALSE, row.names=FALSE, append=FALSE)




