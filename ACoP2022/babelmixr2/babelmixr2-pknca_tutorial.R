#######################################
#
# Installation of required packages
#
#######################################

# For Windows, you must install Rtools, first
# install.packages("installr")
# installr::install.Rtools()

# For Mac, you must install the Mac compilers using the instructions here:
# https://github.com/nlmixr2/nlmixr2#mac-compilation-tools-setup

# For everyone, you must have nlmixr2, nlmixr2lib, babelmixr2, and PKNCA
# installed for this tutorial.
# install.packages(c("nlmixr2", "nlmixr2lib", "babelmixr2", "PKNCA", "units", "pmxTools"))

#######################################
#
# Load packages required for the tutorial
#
#######################################

library(nlmixr2)
library(nlmixr2lib)
library(babelmixr2)

#######################################
#
# Setup data for analysis
#
#######################################

dMod <- nlmixr2data::theo_sd

dMod_nozero <- dMod[(dMod$DV != 0 & dMod$EVID == 0) | (dMod$EVID == 101), ]

#######################################
#
# Simple example of PKNCA use: Automatically estimate starting parameters
#
#######################################

# Prepare a basic model:
# * Load the basic model,
# * add inter-individual variability on ka, cl, and vc, and
# * change the model to use additive residual error, only.
prepModelTheo <-
  nlmixr2lib::readModelDb("PK_1cmt_des") %>%
  addEta(c("ka", "cl", "vc")) %>%
  addResErr("addSd")

# Look at the new model function

prepModelTheo

#######################################
#
# Generate new initial estimates for the model
#
#######################################

# Use `est = "pknca"` to choose to set the new initial estimates wtih PKNCA by
# calculating NCA parameters
#
# Use pkncaControl() to set options for calculation
# * use the dataset with initial concentrations of zero included
#   (`ncaData = dMod`)
#      * `ncaData` can also be used to provide data that PKNCA can use, if it is
#        different from the modeled data.  This may be useful if you only want
#        to include a subset of the studies (e.g. ones with dense PK
#        measurement) in the NCA estimation.
# * provide units for the estimation to modify the dependent variable to
#   correctly assign units
#   (`concu = "mg/L", doseu = "mg/kg", timeu = "hr", volumeu = "L/kg"`)

modPknca1cmt <-
  nlmixr2(
    prepModelTheo,
    data = dMod_nozero, est = "pknca",
    control = pkncaControl(ncaData = dMod, concu = "mg/L", doseu = "mg/kg", timeu = "hr", volumeu = "L/kg")
  )

#######################################
#
# Look at the new model and the NCA estimates
#
#######################################

print(modPknca1cmt)

#######################################
#
# Estimate from that model using typical 'nlmixr2' methods.
#
#######################################

modFocei1cmt <- nlmixr2(modPknca1cmt, data = dMod_nozero, est = "focei")

#######################################
#
# Unit conversion is also possible
#
#######################################

# convert the dataset DV units from mg/L to ng/mL.

dMod_ngmL <- dMod
dMod_ngmL$DV <- dMod_ngmL$DV * 1e3
dMod_ngmL_nozero <- dMod_nozero
dMod_ngmL_nozero$DV <- dMod_ngmL_nozero$DV * 1e3

#######################################
#
# Prepare the model
# * Load the basic model,
# * add inter-individual variability on ka, cl, and vc,
# * change the model to use additive residual error, only, and
# * update the additive residual error to ng/mL units (`cpaddSd <- c(0, 1000)`).
#
#######################################

prepModelTheo_ngmL <-
  nlmixr2lib::readModelDb("PK_1cmt_des") %>%
  addEta(c("ka", "cl", "vc")) %>%
  addResErr("addSd") %>%
  ini(cpaddSd <- c(0, 1000))

modPknca1cmt_ngmL <-
  nlmixr2(
    prepModelTheo_ngmL,
    data = dMod_ngmL_nozero, est = "pknca",
    control = pkncaControl(ncaData = dMod_ngmL, concu = "ng/mL", doseu = "mg/kg", timeu = "hr", volumeu = "L/kg")
  )

print(modPknca1cmt_ngmL)

#######################################
#
# Estimate from that model using typical 'nlmixr2' methods.
#
#######################################

modFocei1cmt_ngmL <- nlmixr2(modPknca1cmt_ngmL, data = dMod_ngmL_nozero, est = "focei")
print(modFocei1cmt_ngmL)

#######################################
#
# Pre-calculate NCA with PKNCA
#
# If automatic calculation of NCA with PKNCA does not work for any reason, you
# may also pre-calculate the parameters.  Then, you can provide the PKNCA
# results objects instead of using the automated calculation.
#
#######################################

library(PKNCA)

d_conc <- dMod_ngmL[dMod_ngmL$EVID == 0, ]
d_dose <- dMod_ngmL[dMod_ngmL$EVID == 101, ]
o_conc <- PKNCAconc(d_conc, DV~TIME|ID)
o_dose <- PKNCAdose(d_dose, AMT~TIME|ID)
o_data <- PKNCAdata(o_conc, o_dose)
o_nca <- pk.nca(o_data)

prepModelTheo_ngmL <-
  nlmixr2lib::readModelDb("PK_1cmt_des") %>%
  addEta(c("ka", "cl", "vc")) %>%
  addResErr("addSd") %>%
  ini(cpaddSd <- c(0, 1000))

modPknca1cmt_manual <-
  nlmixr2(
    prepModelTheo_ngmL,
    data = dMod_ngmL_nozero, est = "pknca",
    control = pkncaControl(ncaResults = o_nca, concu = "ng/mL", doseu = "mg/kg", timeu = "hr", volumeu = "L/kg")
  )

###################
# You must calculate the required parameters of cmax.dn, cl.last, and tmax for
# the conversion to work
###################

d_intervals <- data.frame(start=0, end=Inf, tmax=TRUE, cmax.dn=TRUE, cl.last=TRUE)
o_data <- PKNCAdata(o_conc, o_dose, intervals = d_intervals)
o_nca <- pk.nca(o_data)

modPknca1cmt_manual <-
  nlmixr2(
    prepModelTheo_ngmL,
    data = dMod_ngmL_nozero, est = "pknca",
    control = pkncaControl(ncaResults = o_nca, concu = "ng/mL", doseu = "mg/kg", timeu = "hr", volumeu = "L/kg")
  )

# And then, you can run the model with nlmixr2, like normal.

modFocei1cmt_manual <- nlmixr2(modPknca1cmt_manual, data = dMod_ngmL_nozero, est = "focei")

#######################################
#
# Generate a two-compartment PK model
#
# The link between PKNCA and nlmixr2 works for 1-, 2-, and 3-compartment models.
#
#######################################

prepModelTheo2Cmt <-
  nlmixr2lib::readModelDb("PK_2cmt_des") %>%
  addEta(c("ka", "cl", "vc")) %>%
  addResErr("addSd")

# Generate initial conditions for a 2-compartment model using PKNCA

modPknca2cmt <-
  nlmixr2(
    prepModelTheo2Cmt,
    data = dMod_nozero, est = "pknca",
    control = pkncaControl(ncaData = dMod, concu = "mg/L", doseu = "mg/kg", timeu = "hr", volumeu = "L/kg")
  )

# Look at the new model and the NCA estimates

print(modPknca2cmt)

# Estimate from that model using typical 'nlmixr2' methods.

modFocei2cmt <- nlmixr2(modPknca2cmt, data = dMod_nozero, est = "focei")

#######################################
#
# With the support for model generation, any method may be used for automatic
# (or manual) selection of the best models.
#
#######################################

# Compare and choose the best model by AIC:

AIC(
  modFocei1cmt,
  modFocei2cmt
)

#######################################
#
# The babelmixr2 library adds the ability to use NONMEM or Monolix for the
# estimation.  You must have a license and the ability to run them from the
# computer with the R session.  (There are ways that you can also run them on
# another computer, but those are out of scope for the course.)
#
# Optional: compare these outputs to NONMEM (You must have NONMEM setup so that
# you can run it from the current system)
#
#######################################

## Run both models in NONMEM

## To run in NONMEM, use `est = "nonmem"` and setup the
## `nonmemControl(runCommand)` to the correct NONMEM executable for your system.

modNONMEM1cmt <-
  nlmixr2(
    modPknca1cmt, data = dMod_nozero, est = "nonmem",
    control = nonmemControl(runCommand = "c:/nm75g64/run/nmfe75")
  )

modNONMEM1cmt_ngmL <-
  nlmixr2(
    modPknca1cmt_ngmL, data = dMod_ngmL_nozero, est = "nonmem",
    control = nonmemControl(runCommand = "c:/nm75g64/run/nmfe75")
  )

# Compare the model estimates (they are similar)

fixef(modFocei1cmt)
fixef(modNONMEM1cmt)
fixef(modFocei1cmt_ngmL)
fixef(modNONMEM1cmt_ngmL)

# Compare the AIC (they are also similar)

AIC(
  modFocei1cmt,
  modNONMEM1cmt,
  modFocei1cmt_ngmL,
  modNONMEM1cmt_ngmL
)
