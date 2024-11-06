# Install required packages (only once) ----

installFromCran <- FALSE
installFromGitHub <- FALSE

if (installFromCran) {
  install.packages("nlmixr2", dependencies = TRUE)
  install.packages(c("nonmem2rx", "monolix2rx", "babelmixr2"))
} else if (installFromGitHub) {
  remotes::install_github(
    paste0(
      "nlmixr2/",
      c("rxode2", "nlmixr2est", "nlmixr2", "monolix2rx", "nonmem2rx", "babelmixr2")
    ),
    force = TRUE
  )
}
# Load and examine a typical NONMEM population pharmacokinetic model ----

## Load required packages ----

library(nlmixr2)
library(rxode2)
library(nonmem2rx)
library(babelmixr2)

## View a NONMEM model ----

# This just gets a file that you can use. Use your own file, if you want.
ctlFileName <- system.file("mods/cpt/runODE032.ctl", package="nonmem2rx")
resFileName <- system.file("mods/cpt/runODE032.res", package="nonmem2rx")
resFileName

readLines(ctlFileName, warn = FALSE)
head(readLines(resFileName), 70)

## Load a NONMEM model into rxode2 format ----

modOrigNames <- nonmem2rx(resFileName, save=FALSE, determineError=TRUE)

### Look at the loaded model ----

modOrigNames

# See that the model works comparing NONMEM vs rxode2 (the same works for
# monolix2rx)

plot(modOrigNames)

## You can also plot the differences between model prediction by
## subject (on log scale) ----
plot(modOrigNames, page=1, log="y") # page=TRUE shows all subjects

## Or, you can load the model translating the parameter names to
## easier-to-read names when loading. ----

modCleanNames <-
  nonmem2rx(
    resFileName, save = FALSE,
    thetaNames = c("lcl", "lvc", "lq", "lvp", "prop.sd"),
    etaNames = c("eta.cl", "eta.vc", "eta.q","eta.vp"),
    cmtNames = c("central", "perip")
  )

modCleanNames

## Or, you can clean the model compartment and parameter names manually. ----

modUpdatedNames <-
  modOrigNames |>
  rxRename(
    central = CENTRAL,
    perip = PERI
  )

modUpdatedNames

# Convert NONMEM fit to nlmixr2 fit object (does not re-estimate)

nlmixrFitMod <- as.nlmixr2(modUpdatedNames)

nlmixrFitMod


# Sometimes models are not read with the translated residual error
# (NONMEM only, not monolix2rx).  This can be situation can be
# produced manually by using the option `determineError=FALSE`:

modCleanNamesNoError <-
  nonmem2rx(
    resFileName, save = FALSE,
    thetaNames = c("lcl", "lvc", "lq", "lvp", "prop.sd"),
    etaNames = c("eta.cl", "eta.vc", "eta.q","eta.vp"),
    cmtNames = c("central", "perip"),
    determineError=FALSE)

# In this case you might want to copy the model and manually modify
# the residual error.  You would get the original model, which we are
# changing into a function to simulate the process:
f <- as.function(modCleanNames)

# Once you make the changes you can then convert back by using `as.nonmem2rx()`
f <- as.nonmem2rx(f, modCleanNamesNoError)

# This re-runs the model validation so you know that even though you
# changed the model you can convert this to a nlmixr2 fit with
# confidence that you expressing the model correctly.


# Not that once this is in an appropriate format you can use the data to
# re-estimate (or estimate using new data) in `nlmixr2` if you wish:
# (unlike as.nlmixr2 this performs an estimation)
fit <- nlmixr2(f, modCleanNamesNoError$nonmemData, "focei")
