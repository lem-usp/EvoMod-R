# Fixture: a short run of the C simulation (../c-gsl in the evomod repo) with
# N_e = 60, 4 traits, 10 loci and Delta theta = 0.05 0.05 0 0.
out.dir = test_path("fixtures", "output")
run.dir = file.path(out.dir, "c_run")

read_run = function() suppressMessages(capture.output(
  pop <- ReadFolder("c_run", sel.type = "corridor", output.dir = out.dir)))

test_that("ReadParameters parses pop.parameters.txt", {
  par = ReadParameters(run.dir)
  expect_equal(par$n.e, 60)
  expect_equal(par$n.traits, 4)
  expect_equal(par$n.loci, 10)
  expect_equal(par$delta.theta, c(0.05, 0.05, 0, 0))
  expect_false(is.na(par$seed))
})

test_that("ReadFolder reads every generation with traits from the parameters file", {
  capture.output(pop <- ReadFolder("c_run", sel.type = "corridor", output.dir = out.dir))
  n.gen = length(pop$generation)
  expect_equal(n.gen, 10)
  expect_equal(pop$generation, 1:10)
  expect_equal(dim(pop$g.cor[[1]]), c(4, 4))
  expect_equal(unname(diag(pop$g.cor[[n.gen]])), rep(1, 4))
  expect_true(isSymmetric(pop$g.cor[[n.gen]]))
  expect_length(pop$phenotype[[n.gen]], 4)
  expect_equal(pop$selection.strength, 0.05)
  expect_equal(pop$n.e, 60)
})

test_that("covariances rebuilt from the output files match G of the final population", {
  capture.output(pop <- ReadFolder("c_run", sel.type = "corridor", output.dir = out.dir))
  raw = ReadRawPop(run.dir)
  n.gen = length(pop$g.cov)
  expect_equal(raw$generation, pop$generation[n.gen])
  expect_equal(unname(pop$g.cov[[n.gen]]), unname(raw$G), tolerance = 1e-4)
  expect_equal(unname(diag(pop$g.cov[[n.gen]])), pop$g.var[[n.gen]], tolerance = 1e-6)
})

test_that("ReadRawPop decodes y and B", {
  raw = ReadRawPop(run.dir)
  expect_equal(dim(raw$y), c(20, 60))
  expect_equal(dim(raw$B), c(20, 4, 60))
  expect_true(all(raw$B %in% c(0, 1)))
  expect_equal(dimnames(raw$B)[[2]], paste0("trait", 1:4))
})

test_that("statistics and plots run on simulation output", {
  capture.output(pops <- ReadPattern("c_run", sel.type = "corridor", output.dir = out.dir))
  pop = pops[[1]]
  quiet = function(x) suppressWarnings(suppressMessages(capture.output(r <- x))) 

  expect_length(CalcIsoStat(pop$g.cov, Autonomy), 10)
  expect_length(CalcMeanStat(pop$g.cov[1:2], Autonomy, nsk = 50), 2)

  avg = suppressMessages(AVGRatio(pop$p.cor, 0.05, last.gen = FALSE))
  expect_equal(nrow(avg), 3 * 10)
  expect_equal(unique(avg$generation), 1:10)

  omega = diag(4)
  expect_length(CalcCorrOmegaRS(pop$p.cov[1:2], omega), 2)

  for (p in list(AVGRatioPlot(pops, modules = TRUE),
                 AVGRatioPlot(pops, modules = FALSE),
                 PlotCorrs(pop$p.cor),
                 LastGenMultiStatMultiPlot(pops, function(x) EigenVar(x, 4), "eig"),
                 LastGenStatMultiPlotWithMean(pops, Autonomy, "auto")))
    expect_s3_class(p, "ggplot")
})
