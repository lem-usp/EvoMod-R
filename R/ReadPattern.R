#'Read Pattern from directory
#'
#'Read many populations from a folder using a pattern
#'@param pattern regular expression matched against folder names in \code{output.dir}
#'@param n.traits number of traits; read from each pop.parameters.txt if NULL
#'@param sel.type label stored in \code{selection.type}
#'@param direct.sel passed to \code{ReadFolder}
#'@param output.dir directory holding the simulation output folders
#'@export
ReadPattern <- function(pattern,
                        n.traits = NULL,
                        sel.type,
                        direct.sel = T,
                        output.dir = "./output"){
  folders  <- dir(output.dir, pattern)
  main.data = llply(folders, function(x) ReadFolder(x, n.traits, sel.type, direct.sel, output.dir))
  names(main.data) = folders
  return(main.data)
}
