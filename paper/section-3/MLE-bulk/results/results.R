
#obtain files
type = "ar_0.5";
B = (1:120)

#obtain all the coverage proportions
cov_t = NULL; cov_r = NULL;

dir <- "~/Desktop/Topics/Finished projects/LogisticCov/Code/paper/"
for(b in B){
	newfile = paste0(dir, 'section-3/MLE-bulk/output/',type, '/pct',b,'.txt')
	newcov = as.matrix(read.table(newfile, header=FALSE))
	
	cov_t = rbind(cov_t, newcov)
}

cov_t = cov_t[1:100000,]

#compute statistics
apply(cov_t, 2, mean)*100
apply(cov_t, 2, sd)/sqrt(100000)*100



