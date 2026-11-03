## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
has_desolve <- requireNamespace("deSolve", quietly = TRUE)

## ----library------------------------------------------------------------------
library(r2sbml)
ex <- function(f) system.file("examples", f, package = "r2sbml")

# convertReactions() reports progress on the console; keep only the warnings
convert <- function(file, format = "R", out = tempfile()) {
  tryCatch({
    invisible(capture.output(convertReactions(file, out, format = format)))
    cat("no warning\n")
  }, warning = function(w) cat(conditionMessage(w), "\n"))
  invisible(out)
}
sbml_file <- function(body) {
  f <- tempfile(fileext = ".xml")
  writeLines(sprintf('<?xml version="1.0" encoding="UTF-8" ?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
<model>%s</model>
</sbml>', body), f)
  f
}

## ----algebraic-linear---------------------------------------------------------
linear <- sbml_file('
<listOfParameters><parameter id="k" constant="false"/></listOfParameters>
<listOfRules><algebraicRule><math xmlns="http://www.w3.org/1998/Math/MathML">
<apply><minus/><ci>k</ci><cn>0.9</cn></apply></math></algebraicRule></listOfRules>')
out <- convert(linear, "mrgsolve")
grep("k = ", readLines(out), value = TRUE)

## ----dae-works, eval = has_desolve--------------------------------------------
out <- convert(ex("sbmlalgebraicrules.xml"))
env <- new.env(); source(out, local = env)
d <- as.data.frame(deSolve::daspk(y = env$InitialAmounts, dy = env$InitialDerivatives,
                                  times = seq(0, 10, by = 2), res = env$massBalances,
                                  parms = env$parameters))
max(abs(d$E + d$ES - d$E_total))    # constraint residual, ~1e-13

## ----dae-warns----------------------------------------------------------------
convert(ex("sbmlalgebraicrules.xml"), "mrgsolve")

## ----volume-rule--------------------------------------------------------------
volume <- function(rule) sbml_file(sprintf('
<listOfCompartments><compartment id="c" size="1" constant="false"/></listOfCompartments>
<listOfSpecies><species id="A" compartment="c" initialConcentration="2"
  hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/></listOfSpecies>
<listOfParameters><parameter id="g" value="1" constant="true"/></listOfParameters>
<listOfRules><assignmentRule variable="c">
<math xmlns="http://www.w3.org/1998/Math/MathML">%s</math></assignmentRule></listOfRules>', rule))
time <- '<csymbol encoding="text" definitionURL="http://www.sbml.org/sbml/symbols/time">time</csymbol>'

growing <- volume(sprintf('<apply><plus/><cn>1</cn><apply><times/><ci>g</ci>%s</apply></apply>', time))
out <- convert(growing)
grep("dA_dt = ", readLines(out), value = TRUE)    # V = 1 + g t, so dV/dt = g

## ----volume-rule-warns--------------------------------------------------------
capped <- volume(sprintf('<apply><min/><cn>5</cn><apply><plus/><cn>1</cn>
  <apply><times/><ci>g</ci>%s</apply></apply></apply>', time))
convert(capped)

## ----delay--------------------------------------------------------------------
out <- convert(ex("sbmldelay.xml"))
grep("delay", readLines(out), value = TRUE)

## ----delay-ubiquity-----------------------------------------------------------
convert(ex("sbmldelay.xml"), "ubiquity")

## ----events-------------------------------------------------------------------
out <- convert(ex("sbmlevent.xml"))
grep("events", readLines(out), value = TRUE)

## ----events-inspect-----------------------------------------------------------
getEventMath(getModel(ex("sbmlevent.xml")))

## ----maths--------------------------------------------------------------------
maths <- sbml_file('
<listOfParameters>
<parameter id="x" value="8" constant="true"/>
<parameter id="a" constant="false"/><parameter id="b" constant="false"/>
</listOfParameters>
<listOfRules>
<assignmentRule variable="a"><math xmlns="http://www.w3.org/1998/Math/MathML">
<apply><log/><logbase><cn>2</cn></logbase><ci>x</ci></apply></math></assignmentRule>
<assignmentRule variable="b"><math xmlns="http://www.w3.org/1998/Math/MathML">
<apply><root/><degree><cn>3</cn></degree><ci>x</ci></apply></math></assignmentRule>
</listOfRules>')
out <- convert(maths)
grep("^ +[ab] = ", readLines(out), value = TRUE)
out <- convert(maths, "mrgsolve")
grep("^double [ab] = ", readLines(out), value = TRUE)

## ----undefined----------------------------------------------------------------
undefined <- sbml_file('
<listOfParameters><parameter id="k" constant="true"/></listOfParameters>')
convert(undefined)

## ----reserved-----------------------------------------------------------------
reserved <- sbml_file('
<listOfParameters><parameter id="Inf" value="2" constant="true"/>
<parameter id="v" constant="false"/></listOfParameters>
<listOfRules><assignmentRule variable="v"><math xmlns="http://www.w3.org/1998/Math/MathML">
<apply><times/><ci>Inf</ci><cn>3</cn></apply></math></assignmentRule></listOfRules>')
out <- tempfile()
msg <- capture.output(convertReactions(reserved, out, format = "Julia"))
grep("Renamed", msg, value = TRUE)

## ----species-ic---------------------------------------------------------------
model <- getModel(ex("sbmlsimple.xml"))
getSpeciesIC(model)

## ----species-ic-vs-converted--------------------------------------------------
out <- convert(ex("sbmlsimple.xml"))
grep("^ +E = [0-9]", readLines(out), value = TRUE)   # 5e-07 = 5e-21 / 1e-14

## ----species-table------------------------------------------------------------
getSpeciesTable(model)[, c("ID", "InitialConcentration", "InitialAmount")]

## ----species-names------------------------------------------------------------
getSpeciesNames(model)
getCmtNames(model)

## ----stopping, error = TRUE---------------------------------------------------
try({
getRuleMath(model)     # sbmlsimple.xml has no rules
})

