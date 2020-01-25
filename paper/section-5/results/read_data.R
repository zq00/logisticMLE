dir <- "~/Desktop/Topics/Finished projects/LogisticCov/Code/paper/"
fileloc = "section-5/output/"
names = c("param", "beta_est", "cov_est", "std_est", "dev","rho")

for(name in names){
  obj <- NULL
  
  file_names <- list.files(paste0(dir, fileloc), name)
  for(file in file_names){
    new_name = paste0(dir, fileloc, file)
    
    new_obj <- as.matrix(read.table(new_name,quote="\"", comment.char=""))
    obj <- rbind(obj, new_obj)
    
    assign(name, get("obj")) 
  }
}

tau_hat_1 <- cov_est
tau_hat_2 <- as.matrix(read.table(paste0(dir, fileloc, "hat-tau-rho_1.txt"), quote="\"", comment.char=""))

