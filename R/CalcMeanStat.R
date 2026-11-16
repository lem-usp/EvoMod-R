#'Calculate mean statistic
#'
#'Calculates the mean of a statistic over random directions given a matrix
#'list and an evolqg statistic such as \code{Autonomy}, called as
#'\code{Stat(cov.matrix, beta.mat)}.
#'@export
CalcMeanStat  <- function(mat.list, Stat, nsk = 1000){
  n.traits = dim(mat.list[[1]])[1]
  beta.mat <- array (rnorm (n.traits * nsk), c(n.traits, nsk))
  beta.mat <- apply (beta.mat, 2, Normalize)
  out = lapply(mat.list, function(mat) mean(Stat(cov.matrix = mat, beta.mat = beta.mat)))
  return(unlist(out))
}
