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

## Or, you can load the model translating the parameter names to easier-to-read names when loading. ----

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
