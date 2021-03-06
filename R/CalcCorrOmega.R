#'Caclulate Omega Correlation
#'
#'Calculates the fitness surface correlation with a list of matrices.
#'@export
CalcCorrOmega <- function(mat.list, omega){
  n.traits = dim(mat.list[[1]])[1]
  omega = as.matrix(read.table ("input/omega.csv", header=F, sep=' '))[1:n.traits, 1:n.traits]
  omega = omega[upper.tri(omega)]
  corr.omega <- lapply(mat.list, function(x) cor(x[upper.tri(x)], omega))
  return(unlist(corr.omega))
}

#'Caclulate Omega Correlation
#'
#'Calculates the fitness surface correlation with a list of matrices.
#'@export
CalcCorrOmegaKrz <- function(mat.list, omega){
  n.traits = dim(mat.list[[1]])[1]
  omega = as.matrix(read.table ("input/omega.csv", header=F, sep=' '))[1:n.traits, 1:n.traits]
  corr.omega <- lapply(mat.list, function(x) KrzCor(x, omega, ret.dim = 2))
  return(unlist(corr.omega))
}

#'Caclulate Omega Correlation
#'
#'Calculates the fitness surface correlation with a list of matrices.
#'@export
CalcCorrOmegaRS <- function(mat.list, omega){
  n.traits = dim(mat.list[[1]])[1]
  omega = as.matrix(read.table ("input/omega.csv", header=F, sep=' '))[1:n.traits, 1:n.traits]
  corr.omega <- lapply(mat.list, function(x) RandomSkewers(x, omega)[1])
  return(unlist(corr.omega))
}

#'Caclulate Omega Correlation
#'
#'Calculates the fitness surface correlation with a list of matrices.
#'@export
CalcCorrOmegaEigenVector <- function(mat.list){
  n.traits = dim(mat.list[[1]])[1]
  omega = as.matrix(read.table ("input/omega.csv", header=F, sep=' '))[1:n.traits, 1:n.traits]
  eVec_omega = eigen(omega)$vectors
  corr.omega <- ldply(mat.list, function(x) {
    eVecs = eigen(x)$vectors
    data.frame(PC1 = abs(eVecs[,1]%*%eVec_omega[,1]),
               PC2 = abs(eVecs[,2]%*%eVec_omega[,2]))
  })
  return(corr.omega)
}

#'Caclulate Omega Correlation
#'
#'Calculates the fitness surface correlation with a list of matrices.
#'@export
EigenVar <- function(mat.list, n_eigen){
  n.traits = dim(mat.list[[1]])[1]
  corr.omega <- ldply(mat.list, function(x) eigen(x)$value[1:n_eigen]/sum(diag(x)))
  return(corr.omega)
}


