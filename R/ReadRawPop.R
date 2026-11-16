#' Read original data from dat files
#' 
#' Create R objects of the simulated populations for a single generation
#' from the pop.pop file written by the C simulation.
#' @param folder simulation output folder containing pop.pop and pop.parameters.txt
#' @param label optional label stored in the returned object
#' @import Matrix
#' @export
ReadRawPop <- function(folder, label = NULL){
  parameters = ReadParameters(folder)
  Ne = parameters$n.e
  m = parameters$n.loci
  p = parameters$n.traits
  # pop.pop holds, per individual, y (2m values) then B (p x 2m, row major),
  # followed by the generation number at the very end of the file.
  raw = scan(file.path(folder, "pop.pop"), numeric(), quiet = TRUE)
  rawPop <- array(raw[1:(2*m*(p+1)*Ne)], dim = c(2*m, p+1, Ne))
  dimnames(rawPop) <- list(paste0(c("loci_1-", "loci_2-"), rep(1:m, each = 2)),
                        c("y", paste0("trait", 1:p)),
                        paste0("ind", 1:Ne))
  x <- t(apply(rawPop, 3, function(x) t(x[,-1]) %*% x[,1]))          
  pop = list(y = rawPop[,1,], 
             B = rawPop[,-1,], 
             x = x, 
             G = cov(x), 
             folder = folder, 
             Ne = Ne, 
             m = m, 
             p = p,
             generation = raw[length(raw)],
             label = label)
  
  class(pop) = "sim_pop"
  return(pop)
}

#' @exportS3Method base::print
print.sim_pop <- function(x, ...){
  print(paste("Population Size =", x$Ne))
  print(paste("Number of Loci =", x$m))
  print(paste("Number of traits =", x$p))
  print(paste("Folder =", x$folder))
  if(!is.null(x$label)) print(paste("Label =", x$label)) 
}
