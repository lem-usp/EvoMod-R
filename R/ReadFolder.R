#'Read simulation parameters
#'
#'Read the pop.parameters.txt file written by the C simulation and return
#'a named list with the numeric parameters.
#'@param input.folder simulation output folder
#'@export
ReadParameters <- function(input.folder){
  parameters = scan(file.path(input.folder, "pop.parameters.txt"), character(), quiet = TRUE)
  GetPar = function(name, offset = 2){
    index = which(name == parameters)
    if(length(index) == 0) return(NA)
    as.numeric(parameters[index[1] + offset])
  }
  n.traits = GetPar("n.traits")
  delta.theta = NA
  if("theta" %in% parameters){
    start = which("theta" == parameters)[1] + 2
    delta.theta = as.numeric(parameters[start:(start + n.traits - 1)])
  }
  list(n.e = GetPar("N_e"),
       n.traits = n.traits,
       n.loci = GetPar("n.loci"),
       mu = GetPar("mu"),
       mu.b = GetPar("mu_B"),
       sigma = GetPar("sigma"),
       v.e = GetPar("v_e"),
       burn.in = GetPar("burn_out"),
       stabilizing = GetPar("stabilizing"),
       selective = GetPar("selective"),
       delta.theta = delta.theta,
       seed = GetPar("seed"))
}

#'Read Folder
#'
#'Read a single population folder and return a list with all poplation atributes
#'@param input.folder name of the folder inside \code{output.dir}
#'@param n.traits number of traits; read from pop.parameters.txt if NULL
#'@param sel.type label stored in \code{selection.type}
#'@param direct.sel if TRUE, \code{selection.strength} is the first entry of Delta theta, otherwise 0
#'@param output.dir directory holding the simulation output folders
#'@export
ReadFolder  <- function(input.folder, n.traits = NULL, sel.type, direct.sel = T,
                        output.dir = "./output"){
  print(input.folder)
  input.folder = file.path(output.dir, input.folder)
  parameters = ReadParameters(input.folder)
  if(is.null(n.traits))
    n.traits = parameters$n.traits
  p.cor = ReadMatrices(file.path(input.folder, "p.corr.dat"), n.traits)
  g.cor = ReadMatrices(file.path(input.folder, "g.corr.dat"), n.traits)
  p.var = ReadVariances(file.path(input.folder, "p.var.dat"), n.traits)
  g.var = ReadVariances(file.path(input.folder, "g.var.dat"), n.traits)
  h.var = ReadVariances(file.path(input.folder, "h.var.dat"), n.traits)
  phenotype = ReadVariances(file.path(input.folder, "phenotype.dat"), n.traits)
  input.file = file.path(input.folder, "b.mat.dat")
  if(file.exists(input.file))
    b.mat = ReadBMat(input.file, n.traits, n.loci = 2*parameters$n.loci)
  else
    b.mat = NULL

  p.cov = CalcCovar(p.cor, p.var)
  g.cov = CalcCovar(g.cor, g.var)

  if(direct.sel)
    selection.strength = parameters$delta.theta[1]
  else
    selection.strength = 0.

  out.list = list(p.cor = p.cor,
                  g.cor = g.cor,
                  p.var = p.var,
                  g.var = g.var,
                  h.var = h.var,
                  p.cov = p.cov,
                  g.cov = g.cov,
                  phenotype = phenotype,
                  b.mat = b.mat,
                  mu.b = parameters$mu.b,
                  n.e = parameters$n.e,
                  parameters = parameters,
                  selection.type = sel.type,
                  selection.strength = selection.strength,
                  generation = as.numeric(names(p.var)))
  return(out.list)
}
