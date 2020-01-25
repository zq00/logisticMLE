dir <- "~/Desktop/Topics/Finished projects/LogisticCov/Code/paper/"
fileloc = "genetics/output/"
names = c("param", "beta_est", "cov_est", "std_est", "dev","gamma_hat","pc")

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

