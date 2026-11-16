#'AVG Ratio Calculation
#'
#'Calculate the AVG ratio for a list of matrices, returning a data.frame
#'@param mat.list named list of correlation matrices (names are generations)
#'@param Selection_Strength value stored in the Selection_Strength column
#'@param last.gen if TRUE, only the last matrix is used
#'@param generations generation of each matrix; defaults to \code{names(mat.list)}
#'@export
AVGRatio <- function(mat.list, Selection_Strength, last.gen = T, generations = NULL){
  n.traits = dim(mat.list[[1]])[1]
  module.1 = rep(c(1,0), each = n.traits/2)
  module.2 =  rep(c(0,1), each = n.traits/2)
  modularity.hipot = cbind(module.1, module.2)
  if(last.gen)
    AVGRatio <- ldply(mat.list[length(mat.list)], function(x) TestModularity(x, modularity.hipot))
  else{
    if(is.null(generations))
      generations = as.numeric(names(mat.list))
    if(length(generations) != length(mat.list))
      stop("generations must have one entry per matrix")
    AVGRatio <- ldply(mat.list, function(x) TestModularity(x, modularity.hipot))
    n.hipot = nrow(AVGRatio) / length(mat.list)
    AVGRatio['generation'] = rep(generations, each = n.hipot)
  }
  AVGRatio['Selection_Strength'] = Selection_Strength
  return(AVGRatio)
}

#'Simple AVG Ratio Calculation
#'
#'Calculate the AVG ratio for a list of matrices, returning a data.frame
#' @export
AVGRatioSimple <- function(mat.list, num.cores = 4){
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
    parallel = TRUE
  }
  else{
    parallel = FALSE
  }
  n.traits = dim(mat.list[[1]])[1]
  module.1 = rep(c(1,0), each = n.traits/2)
  module.2 =  rep(c(0,1), each = n.traits/2)
  modularity.hipot = cbind(module.1, module.2)
  cor.y = CreateHipotMatrix(modularity.hipot)[[3]]
  index <- cor.y[lower.tri(cor.y)]
  AVGRatioWrap = function(cor.x){
    avg.plus <- mean (cor.x [lower.tri(cor.x)] [index != 0])
    avg.minus <- mean (cor.x [lower.tri(cor.x)] [index == 0])
    avg.ratio <- avg.plus / avg.minus
    return(avg.ratio)
  }
  avg.ratio = laply(mat.list, AVGRatioWrap, .parallel=parallel)
  return(avg.ratio)
}