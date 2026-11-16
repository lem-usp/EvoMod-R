#'Plot Correlations
#'
#'Plot within module and between module average correlations through time
#'@param p.cor named list of correlation matrices (names are generations)
#'@param generations generation of each matrix; defaults to \code{names(p.cor)}
#'@export
PlotCorrs = function(p.cor, generations = as.numeric(names(p.cor))){
  avg = AVGRatio(p.cor, NA, last.gen = FALSE, generations = generations)
  m.avg = melt(avg[, c('hypothesis', 'generation', 'AVG+', 'AVG-')],
               id.vars = c('hypothesis', 'generation'))
  m.avg = m.avg[(m.avg$hypothesis != "Full Integration" & m.avg$variable == "AVG+") |
                (m.avg$hypothesis == "Full Integration" & m.avg$variable == "AVG-"), ]
  labels = c("Within Module 1", "Within Module 2", "Between Modules")
  m.avg$group = factor(paste(m.avg$hypothesis, m.avg$variable),
                       levels = c("module.1 AVG+", "module.2 AVG+", "Full Integration AVG-"),
                       labels = labels)
  avg.plot = ggplot(m.avg, aes(generation, value, colour = group, shape = group)) +
    geom_point(alpha = 0.65) +
    labs(x="Generation",
         y="Average Correlation",
         colour = "Module", shape = "Module") + theme_bw()
  return(avg.plot)
}
