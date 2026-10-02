#################################################################################
# ---- Exercises ----------------------------------------------
##  Hands-on bonus assignment: ODE update                                      ##
##  -Change the system of ODEs and examine the results                         ##
##  For example: extend the original model with five transit compartments      ##
##  and use 4 bolus doses in the 1st compartment                               ##
#################################################################################

library(nlmixr2)     # model simulation (rxSolve, event tables)
library(tidyverse)   # data handling (tibble, mutate, filter, group_by, pivot_longer) and ggplot2
library(xgxr)        # helper functions for exploratory PK/PD plots (not used further in this script)
library(ggPMX)       # diagnostic plots for pharmacometric models (not used further in this script)

# Use all available threads for rxode2 computations
setRxThreads(percent=100)


# Load the saved fit run101 (used for the simulations below)
load("run101.Rdata")


# ---- 1.single subject simulation --------------------------------------------

## ---- 1.1 single dose administration ----------------------------------------
# A 500 mg initial dose into compartment 1


ev <- eventTable() ## create an empty event table
ev$add.dosing(dose = 500) # add a single dose of 500 mg (into compartment 1, the default)

## add time points to the event table where concentrations will be simulated 
## these actions are cumulative 

ev$add.sampling(seq(0, 120, 0.1))  # every 0.1 h from 0 to 120 h

# simulate a single subject with the model run101 (loaded above) and the event table
sim_typ1 <- rxSolve(run101, ev)    

# plot function in rxode2: the amount in the depot compartment and the concentration cp over time
plot(sim_typ1, depot, cp)

# but also possible to use ggplot
# the simulation is returned as a tibble and converted to the long format:
# one row per time and variable (all columns except time), the variable name is stored in CMT
sim_typ_1l <- rxSolve(run101, ev,  returnType = "tibble" )  %>% 
  pivot_longer(cols = !time, values_to = "PRED", names_to = "CMT") 



## then plot the simulated outcomes:
## all simulated variables (compartment amounts and concentration), each in its own panel

# plot all compartments

ggplot(data  =sim_typ_1l, aes(x= time, y  = PRED)) + 
  geom_line() +
  facet_wrap(~CMT)

# or filter what is needed (here only the concentration cp)
ggplot(data  = sim_typ_1l |>
         filter(CMT == "cp" ), aes(x= time, y  = PRED)) + 
  geom_line() 


## Extend the eventTable by adding three infusions to the central compartment
## Remember: updates to the eventTable are cumulative
## Add three 2-hour infusions of 250 mg every 12 h, starting at 36 h,
## into compartment 2, simulated over 120 h.
ev$add.dosing(
  dose = 250,           #mg
  nbr.doses = 3,        #add three doses
  dosing.to = 2,        #add them to the second ODE in the model (=central)
  dosing.interval = 12, #h; set the doses 12 hours apart
  rate = 125,           #mg/h; infuse at a rate of 125 mg/h, resulting in 2-hour infusions
  start.time = 36       #h; have the three doses start at 36h
)


# simulate again with the extended event table (500 mg dose + three infusions)
sim_typ2 <- rxSolve(run101, ev )



## then plot the simulated outcomes:
## the amount in the depot compartment and the concentration cp

plot(sim_typ2, depot, cp)

# ---- 1. use a data frame structure as input ---------------------------------------------------

##  Hands-on bonus assignment: ODE update                                      ##
##  -Change the system of ODEs and examine the results                         ##
##  For example: extend the original model with five transit compartments      ##
##  and use 4 bolus doses in the 1st compartment                               ##
