#'Calculate v-cv matrix
#'
#'Calculate a v-cv matrix using a list of correlations and variances.
#'The simulation writes variances (diagonals of G and P), so covariances are
#'\code{cor * sqrt(var_i * var_j)}.
#'@export
CalcCovar  <- function(corr, vars){
  num.gens = length(vars)
  covs = vector('list', num.gens)
  for(i in 1:num.gens)
    covs[[i]] = corr[[i]]*sqrt(outer(vars[[i]], vars[[i]]))
  names(covs) = names(vars)
  return(covs)
}
