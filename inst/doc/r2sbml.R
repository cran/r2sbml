## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
has_desolve <- requireNamespace("deSolve", quietly = TRUE)

## ----library------------------------------------------------------------------
library(r2sbml)

## ----examples-----------------------------------------------------------------
basename(list.files(system.file("examples", package = "r2sbml"), pattern = "[.]xml$"))

## ----getmodel-----------------------------------------------------------------
sbml_file <- system.file("examples", "sbmlsimple.xml", package = "r2sbml")
model <- getModel(sbml_file)
class(model)

## ----tables-------------------------------------------------------------------
getSpeciesTable(model)

getReactionTable(model)

## ----vectors------------------------------------------------------------------
getCmtNames(model)
getCmtSizes(model)
getNumSpecies(model)

## ----printing-----------------------------------------------------------------
rules_model <- getModel(system.file("examples", "sbmlassignmentrules.xml",
                                    package = "r2sbml"))
getRuleMath(rules_model)

## ----capture------------------------------------------------------------------
txt <- capture.output(getReactionMath(model))
txt

## ----stopping, error = TRUE---------------------------------------------------
try({
# sbmlsimple.xml has no rules at all
getRuleMath(model)
})

## ----looping------------------------------------------------------------------
safely <- function(expr) tryCatch(expr, error = function(e) NULL)
files  <- list.files(system.file("examples", package = "r2sbml"),
                     pattern = "[.]xml$", full.names = TRUE)

has_rules <- vapply(files, function(f) {
  !is.null(safely(capture.output(getRuleMath(getModel(f)))))
}, logical(1))
data.frame(model = basename(files), has_rules = unname(has_rules))

## ----convert-basic------------------------------------------------------------
out_r <- tempfile(fileext = ".R")
convertReactions(sbml_file, out_r, format = "R")

## ----bad-format, error = TRUE-------------------------------------------------
try({
convertReactions(sbml_file, tempfile(), format = "fortran")
})

## ----show-r-------------------------------------------------------------------
cat(readLines(out_r), sep = "\n")

## ----run-r, eval = has_desolve------------------------------------------------
library(deSolve)

env <- new.env()
source(out_r, local = env)

solution <- ode(y     = env$InitialAmounts,
                times = seq(0, 5, by = 1),
                func  = env$massBalances,
                parms = env$parameters)
round(solution, 9)

## ----conservation, eval = has_desolve-----------------------------------------
# free + bound enzyme is conserved, and so is total substrate
d <- as.data.frame(solution)
range(d$E + d$ES)
range(d$S + d$ES + d$P)

## ----mrgsolve-----------------------------------------------------------------
out_mrg <- tempfile(fileext = ".cpp")
convertReactions(sbml_file, out_mrg, format = "mrgsolve")
cat(readLines(out_mrg), sep = "\n")

## ----rxode--------------------------------------------------------------------
out_rx <- tempfile(fileext = ".R")
convertReactions(sbml_file, out_rx, format = "nlmixr2")
cat(readLines(out_rx), sep = "\n")

## ----matlab-------------------------------------------------------------------
out_m <- tempfile(fileext = ".m")
convertReactions(sbml_file, out_m, format = "MATLAB")
cat(readLines(out_m), sep = "\n")

## ----julia--------------------------------------------------------------------
out_jl <- tempfile(fileext = ".jl")
convertReactions(sbml_file, out_jl, format = "Julia")
cat(readLines(out_jl), sep = "\n")

## ----ubiquity-----------------------------------------------------------------
out_ub <- tempfile(fileext = ".txt")
convertReactions(sbml_file, out_ub, format = "ubiquity")
cat(readLines(out_ub), sep = "\n")

## ----assignment---------------------------------------------------------------
assign_file <- system.file("examples", "sbmlassignmentrules.xml", package = "r2sbml")
out_assign  <- tempfile(fileext = ".R")
convertReactions(assign_file, out_assign, format = "R")

lines <- readLines(out_assign)
cat(lines[grep("Assignment Rules", lines):length(lines)], sep = "\n")

## ----algebraic----------------------------------------------------------------
alg_file <- system.file("examples", "sbmlalgebraicrules.xml", package = "r2sbml")
out_alg  <- tempfile(fileext = ".R")
convertReactions(alg_file, out_alg, format = "R")

lines <- readLines(out_alg)
cat(lines[grep("Residuals:", lines)[1]:(grep("return\\(list\\(Residuals", lines))], sep = "\n")

## ----algebraic-solve, eval = has_desolve--------------------------------------
env <- new.env()
source(out_alg, local = env)

sol <- deSolve::daspk(y     = env$InitialAmounts,
                      dy    = env$InitialDerivatives,
                      times = seq(0, 10, by = 2),
                      res   = env$massBalances,
                      parms = env$parameters)
round(sol, 6)

## ----algebraic-check, eval = has_desolve--------------------------------------
d <- as.data.frame(sol)
# both algebraic rules hold along the whole trajectory
max(abs(d$E + d$ES - d$E_total))
max(abs(1 * d$E * d$S - (0.5 + 0.5) * d$ES))

## ----algebraic-matlab---------------------------------------------------------
out_alg_m <- tempfile(fileext = ".m")
convertReactions(alg_file, out_alg_m, format = "MATLAB")
lines <- readLines(out_alg_m)
cat(grep("M = eye|M\\(|MassSingular|ode15s", lines, value = TRUE), sep = "\n")

## ----algebraic-julia----------------------------------------------------------
out_alg_jl <- tempfile(fileext = ".jl")
convertReactions(alg_file, out_alg_jl, format = "Julia")
lines <- readLines(out_alg_jl)
cat(grep("Solve with|^M |^M\\[|mass_matrix|ODEProblem", lines, value = TRUE), sep = "\n")

## ----growing-model------------------------------------------------------------
growing <- tempfile(fileext = ".xml")
writeLines('<?xml version="1.0" encoding="UTF-8" ?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
<model substanceUnits="mole" volumeUnits="litre" timeUnits="second" extentUnits="mole">
<listOfCompartments>
<compartment id="c" size="1" spatialDimensions="3" constant="false"/>
</listOfCompartments>
<listOfSpecies>
<species id="A" compartment="c" initialConcentration="2" boundaryCondition="false"
hasOnlySubstanceUnits="false" constant="false"/>
</listOfSpecies>
<listOfParameters>
<parameter id="g" value="1" constant="true"/>
</listOfParameters>
<listOfRules>
<rateRule variable="c">
<math xmlns="http://www.w3.org/1998/Math/MathML"><ci>g</ci></math>
</rateRule>
</listOfRules>
</model>
</sbml>', growing)

out_grow <- tempfile(fileext = ".R")
convertReactions(growing, out_grow, format = "R")
lines <- readLines(out_grow)
cat(grep("InitialAmounts|^ +[Ac] =|_dt|states\\[", lines, value = TRUE), sep = "\n")

## ----growing-solve, eval = has_desolve----------------------------------------
env <- new.env()
source(out_grow, local = env)

times <- c(0, 1, 3, 7)
d <- as.data.frame(deSolve::ode(y = env$InitialAmounts, times = times,
                                func = env$massBalances, parms = env$parameters,
                                rtol = 1e-10, atol = 1e-10))
data.frame(time = times, A = d$A, exact = 2 / (1 + times), amount = d$A * d$c)

