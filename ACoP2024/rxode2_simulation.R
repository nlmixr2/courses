## Expanding on the nonmem2rx model, we can simulate directly from the
## imported model.

## We will use the case where the error did not actually get imported
## in a nlmixr2-compatible format.  While we cannot estimate from this
## model, we can simulate from it.

resFileName <- system.file("mods/cpt/runODE032.res", package="nonmem2rx")
resFileName

library(babelmixr2)

set.seed(42)
rxode2::rxSetSeed(42)


modCleanNamesNoError <-
  nonmem2rx(
    resFileName, save = FALSE,
    thetaNames = c("lcl", "lvc", "lq", "lvp", "prop.sd"),
    etaNames = c("eta.cl", "eta.vc", "eta.q","eta.vp"),
    cmtNames = c("central", "perip"),
    determineError=FALSE)


# First lets create an event table for the simulation

ev <- et(amt=120000, time=0, addl=9, ii=12) %>%
  et(amt=120000, time=120, addl=5, ii=24)
# You can see the dosing here
ev

# Now lets add some observations
ev <- ev %>%
  et(0, 264, length.out=500)

# Now that you have the event table, and the model, you can simulate

sim <- rxSolve(modCleanNamesNoError, ev)

# No you can see the results

plot(sim, f) # f= the results you want to see from the model (you can change this)

plot(sim, f, log="y") # You can also change the y-axis to log scale

# You could expand this with a couple of subjects
ev3 <- ev %>%
  et(id=1:3)

sim3 <- rxSolve(modCleanNamesNoError, ev3)

plot(sim3, f)

# If you want to see a simulation of a single study you could simulate
# that as well by specifying the number of subjects (say 40)
ev40 <- ev %>%
  et(id=1:40)

sim40 <- rxSolve(modCleanNamesNoError, ev40)

# You can see all the subjects
plot(sim40, f)

# Or you can see the confidence bands with confint:
ci40 <- confint(sim40, parm="f", level=0.90)
print(ci40)
plot(ci40)

# You can even simulate using the uncertanty in your parameters

ev <- et(amt=120000, time=0, addl=9, ii=12) %>%
  et(amt=120000, time=120, addl=5, ii=24)
# You can see the dosing here
ev


# Simulate the 50 trails of the same design:
sim40.50 <- rxSolve(modCleanNamesNoError, ev40, nStud=50)

# You can see that the population parameters change by taking a summary:

summary(sim40.50$params)

# You can see the simulated values for each population parameter with
# the $thetaMat property:
sim40.50$thetaMat

# You can also see the omegas sampled from the inverse wishart
# distribution used for the simulation:
head(sim40.50$omegaList)

# You can also see the sigmas sampled from the inverse wishart
# distribution used for the simulation:
head(sim40.50$sigmaList)


# These come from the imported simulation parameters in the model
# iteslf, that is `thetaMat` representing the covariance matrix of the
# population parameters, and `dfSub` and `dfObs` representing the
# degrees of freedom for the inverse wishart distribution for the
# omega and sigma matrices, respectively.

# You can also use the confint to produce a confidence interval for
# this simulation:

ci40.50 <- confint(sim40.50, parm="f", level=0.9)

# You can see the confidence intervals in the printouts and plots
print(ci40.50)

plot(ci40.50)
