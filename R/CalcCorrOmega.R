# Fitness surface matrix: the one given, or the C simulation's input/omega.csv
GetOmega <- function(omega, n.traits){
  if(is.null(omega))
    omega = as.matrix(read.table("input/omega.csv", header=F, sep=' '))
  as.matrix(omega)[1:n.traits, 1:n.traits]
}

#'Caclulate Omega Correlation
#'
#'Calculates the fitness surface correlation with a list of matrices.
#'@param mat.list list of matrices
#'@param omega fitness surface matrix; if NULL, read from input/omega.csv
#'@export
CalcCorrOmega <- function(mat.list, omega = NULL){
  n.traits = dim(mat.list[[1]])[1]
  omega = GetOmega(omega, n.traits)
  omega = omega[upper.tri(omega)]
  corr.omega <- lapply(mat.list, function(x) cor(x[upper.tri(x)], omega))
  return(unlist(corr.omega))
}

#'Caclulate Omega Correlation
#'
#'Calculates the fitness surface correlation with a list of matrices.
#'@param mat.list list of matrices
#'@param omega fitness surface matrix; if NULL, read from input/omega.csv
#'@export
CalcCorrOmegaKrz <- function(mat.list, omega = NULL){
  n.traits = dim(mat.list[[1]])[1]
  omega = GetOmega(omega, n.traits)
  corr.omega <- lapply(mat.list, function(x) KrzCor(x, omega, ret.dim = 2))
  return(unlist(corr.omega))
}

#'Caclulate Omega Correlation
#'
#'Calculates the fitness surface correlation with a list of matrices.
#'@param mat.list list of matrices
#'@param omega fitness surface matrix; if NULL, read from input/omega.csv
#'@export
CalcCorrOmegaRS <- function(mat.list, omega = NULL){
  n.traits = dim(mat.list[[1]])[1]
  omega = GetOmega(omega, n.traits)
  corr.omega <- lapply(mat.list, function(x) RandomSkewers(x, omega)[1])
  return(unlist(corr.omega))
}

#'Caclulate Omega Correlation
#'
#'Calculates the fitness surface correlation with a list of matrices.
#'@param mat.list list of matrices
#'@param omega fitness surface matrix; if NULL, read from input/omega.csv
#'@export
CalcCorrOmegaEigenVector <- function(mat.list, omega = NULL){
  n.traits = dim(mat.list[[1]])[1]
  omega = GetOmega(omega, n.traits)
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


