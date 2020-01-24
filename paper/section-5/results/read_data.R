fileloc = "/section-5/output/"
ending = c("_1_1.txt", "_1_2.txt", "_2.txt")
names = c("param", "beta_est", "cov_est", "std_est", "dev","rho")
# because size restriction to uploaded data, beta_est, cov_est, std_est are splitted into 10 files
# you need to adjust the endings here to load them properly

for(name in names){
  obj <- NULL
  
  for(end in ending){
    new_name = paste0(fileloc, name, end)
    
    new_obj <- as.matrix(read.table(new_name,quote="\"", comment.char=""))
    obj <- rbind(obj, new_obj)
    
    assign(name, get("obj")) 
  }
}

tau_hat_1 <- cov_est
tau_hat_2 <- as.matrix(read.table(paste0(fileloc, "hat-tau-rho_1.txt"), quote="\"", comment.char=""))

