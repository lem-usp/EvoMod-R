#'Caculate ISometrics statistic
#'
#'Calculates a statistic in the direction of the first module (first half of
#'the traits) given a matrix list and an evolqg statistic such as
#'\code{Autonomy}, called as \code{Stat(cov.matrix, beta.mat)}.
#'@export
CalcIsoStat  <- function(mat.list, Stat){
  betas = rep(c(1, 0), each=dim(mat.list[[1]])[1]/2)
  betas = matrix(Normalize(betas))
  out = laply(mat.list, function(mat) Stat(cov.matrix = mat, beta.mat = betas))
  return(out)
}
