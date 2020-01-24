#--- Job submission arguments (please ignore)
args = commandArgs(trailingOnly=TRUE)
array_ind = as.numeric(args[1])
task_ind <- as.numeric(args[2])
ind = array_ind * 1000 + task_ind; print(ind)

sim_ind = 2; 
fileloc = paste0("/section-5/output/",sim_ind,"/data_",ind,"/") 
print(fileloc)
system(paste0('mkdir -p ',fileloc)); 
# --- End of job specific parameters

# phase transition curve 
phase0 = read.table("/section-5/code/phase0.txt")

#-- number of samples
BB = 1; 	#number of simulations (X,y)  
#--

#-- parameters
n <- 4000
kappa <- 0.2
p_sig <- 0.5
gamma2 <- 5
#--

p <- n * kappa
alpha_s <- 1.4994
sigma_s <- 4.7436
lambda_s <- 3.0269

#-- Parameters ProbeFrontier
B = 50; #num of repitition

#-- store estimated parameters
gamma_store = numeric(BB); 
param = matrix(0, BB , 3);

#-- read parameters
Sigma <- as.matrix(read.table(paste0("/section-5/param/Sigma.txt")))  
inv_Sigma <- 1/sqrt(diag(solve(Sigma)))
R <- chol(Sigma)

beta <- numeric(p) 
non_null <- scan(paste0("/section-5/param/non_null.txt"))
j_null <- scan(paste0("/section-5/param/jnull.txt"))
beta_solve <- uniroot(function(x){ beta[non_null] <- x; return(t(beta)%*%(Sigma%*%beta)/n - gamma2);}, c(0,300), check.conv=TRUE,extendInt = "upX")$root
beta[non_null] <- beta_solve
print(abs( t(beta)%*%(Sigma%*%beta)/n - gamma2)) #check if desired signal strength

#-- ProbeFrontier estimate of gamma2
kappas = seq(kappa, 0.5, by = 0.01); n_kappa = length(kappas);
pi_hats = numeric(n_kappa);
is_sep = numeric(B)

for(bb in 1:BB){
  #samples one observation
  X <- matrix(rnorm(n*p, 0, 1), n, p) %*% R/ sqrt(n); eta <- X %*% beta; mu <- exp(eta) / (1+exp(eta))
  y <- rbinom(n, 1, mu)
  
  #write data		
  write.table(X,paste0(fileloc,"X.txt"),col.names = FALSE, row.names=FALSE);
  write.table(y,paste0(fileloc,"Y.txt"),col.names = FALSE, row.names=FALSE);
  
  #-- Begin ProbeFrontier
  for(i in 1:n_kappa){
    stop=FALSE;
    for(b in 1:B){
      k = kappas[i]; nn = floor(p / k); #number of samples in subsample; changed dec.3
      
      sind = sample(1:n, nn, replace = FALSE); #subsample without replacement
      Xs = X[sind,];Ys = y[sind];
      write.table(Xs,paste0(fileloc,'Xs.txt'),col.names = FALSE, row.names=FALSE);
      write.table(Ys,paste0(fileloc,'Ys.txt'),col.names = FALSE, row.names=FALSE);
      
      #check if mle exists
      command = paste0('ml load matlab; matlab -nodisplay -r "addpath /section-5/code/; fileloc = \'', fileloc, ' \';is_sep; exit;"') # run the matlab code to calculate parameters
      system(command)
      #obtain result
      is_sep[b] = scan(paste0(fileloc,'sep.txt'))
      print(c(b,is_sep[1:b]));
    }
    pi_hats[i] = mean(is_sep);
    print(c(i,pi_hats[1:i]));
    if(pi_hats[i]>=0.5){ 
      stop = TRUE; #put a flag that have found kappa			
      
      kappa_n = seq(kappas[i-1], kappas[i], by=0.001) #check more values
      pi_n = numeric(length(kappa_n))
      
      pi_n[1] = pi_hats[i-1]; pi_n[length(kappa_n)] = pi_hats[i]
      kappa_hat = kappas[i];	
      
      for(j in 2:(length(kappa_n)-1)){
        k = kappa_n[j]; nn = floor(p/k); #changed dec.3
        for(b in 1:B){		
          sind = sample(1:n, nn, replace = FALSE);
          Xs = X[sind,];Ys = y[sind];
          write.table(Xs,paste0(fileloc,'Xs.txt'),col.names = FALSE, row.names=FALSE);
          write.table(Ys,paste0(fileloc,'Ys.txt'),col.names = FALSE, row.names=FALSE);	
          #check if mle exists
          command = paste0('ml load matlab; matlab -nodisplay -r "addpath /section-5/code/; fileloc = \'', fileloc, ' \';is_sep; exit;"')
          system(command)
          #obtain result
          is_sep[b] = scan(paste0(fileloc,'sep.txt'))
          print(c(b,is_sep[1:b]));
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
  gamma_store[bb] = gamma_hat;
  write.table(gamma_hat,paste0(fileloc,'gamma.txt'),col.names=FALSE,row.names=FALSE);
  
  command = paste0('ml load matlab; matlab -nodisplay -r "addpath /section-5/code/; fileloc = \'', fileloc, ' \';kappa=0.2; FindParam; exit;"')
  system(command);
  
  param_new = as.matrix(read.table(paste0(fileloc,"param.txt")));
  param[bb,] = param_new[1,]; 
  
  print(param);
  #-- end of ProbeFrontier
  
  #-- Run Logistic regression on (X,y)
  cov_store <- 1/sqrt(diag(solve(t(X) %*% X)))
  
  #fit a logistic regression
  fit <- glm(y~X+0, family=binomial)
  beta_est <- fit$coefficients
  std_est <- summary(fit)$coef[,2]
  
  fit0 <- glm(y~X[,-j_null]+0, family=binomial)
  dev_store <- -anova(fit, fit0)$Deviance[2]
  
  ci_est <- cbind( fit$coefficients - 1.96 * sigma_s/(cov_store/ sqrt(1-kappa)),fit$coefficients +  1.96 * sigma_s/(cov_store / sqrt(1-kappa)))/alpha_s
  ci_t <-  cbind( fit$coefficients - 1.96 * sigma_s/(inv_Sigma),fit$coefficients +  1.96 * sigma_s/(inv_Sigma))/alpha_s
  ci_r <- cbind( fit$coefficients - 1.96 * std_est,fit$coefficients +  1.96 * std_est)
  
  #compute coverage proportion
  pc_t <- mean((ci_t[,2] >= beta) & (ci_t[,1] <= beta))
  pc_est <- mean((ci_est[,2] >= beta) & (ci_est[,1] <= beta))
  pc_r <- mean((ci_r[,2] >= beta) & (ci_r[,1] <= beta))
  
  std_store = std_est
  print(c(pc_est,pc_t,pc_r))
  #ProbeFrontier estimaties
  alpha_hat = param_new[1,1]
  sigma_hat = param_new[1,2]
  lambda_hat = param_new[1,3]
  
  ci_pf = cbind(fit$coefficients - 1.96 * sigma_hat/(cov_store/(sqrt(1-kappa))), fit$coefficients + 1.96 * sigma_hat / (cov_store/sqrt(1-kappa))) / alpha_hat
  pc_pf = mean((ci_pf[,2]>=beta) & (ci_pf[,1]<=beta))

  print(paste0("simulation",bb,"coverage by ProbeFrontier is", pc_pf));
}

#-- write results
write.table(param,paste0(fileloc,"param_",ind,'.txt'),col.names=FALSE, row.names=FALSE,append=FALSE)

