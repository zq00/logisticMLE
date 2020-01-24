#randomly generate a covariance matrix
#reference: How to generate random matrices from the classical compact groups, Francesco Mezzadri

#--- File location
fileloc = '/section-3/param/'

#--- Parameters
p = 800; #number of variables

create.Sigma <- function(p, type, ...){
  param = list(...)
  if(type=='random'){
    # randomly sample eigenvalues
    eig <- rchisq(p, df = 10)
    # randomly sample a rotation matrix U
    Z = matrix(rnorm(p^2), p)
    #take QR recomposition
    U = qr.Q(qr(Z))
    
    # covariance matrix
    Sigma <- t(U) %*% (diag(eig) %*% U) 
    d <- diag(Sigma)
    # correlation matrix
    Sigma / sqrt(d %*% t(d))
    
  }else if(type == 'ar'){
    #sample a covariance matrix from AR(1) model
    toeplitz(param$rho^(0:(p-1)))
    
  }else if(type == 'identity'){
    diag(rep(1, times = p))
  }
}

subfolder <- "random/"

Sigma <- create.Sigma(p, "random")

write.table(Sigma, paste0(fileloc, subfolder, '/Sigma.txt'), col.names=FALSE, row.names=FALSE)

#sample non-null coordinates
non_null = sample(1:p, p/2, replace=FALSE)
write.table(non_null, paste0(fileloc, subfolder, 'non_null.txt'), col.names=FALSE, row.names=FALSE)

#sample one null
j_nonnull = sample((1:p)[non_null], 10, replace = FALSE)
write.table(j_nonnull, paste0(fileloc, subfolder, 'j_nonnull.txt'), col.names=FALSE, row.names=FALSE)
