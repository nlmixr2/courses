#################################################################################
# ---- Exercises ----------------------------------------------
##  Examine the GOF plots and implement alternative absorption models          ##
##  - one or more transit compartment(s)                                       ##
##                                                                             ##
##  Compare vpcs of alternatives and compare OFVs:                             ##
##                                                                             ##
##  run100$objf - run101$objf                                                  ##
##                                                                             ##
#################################################################################




library(nlmixr2)     # model definition and estimation
library(tidyverse)   # loaded here; not used further in this script
library(xgxr)        # loaded here; not used further in this script
library(ggPMX)       # diagnostic (GOF) plots of the fitted models

# Use all available threads for rxode2 computations
setRxThreads(percent=100)


# Let's read the data 

PKdata <- read.csv("data/warfarin_PKS.csv")


# ---- 1. Base model: one compartment, first-order absorption -------------------

One.comp.KA.ODE <- function() {
  ini({
    # Where initial estimates are specified
    lka  <- log(1.15)  #log ka (1/h)
    lcl  <- log(0.135) #log Cl (L/h)
    lv   <- log(8)     #log V (L)
    add.err  <- 0.6    #additive error (mg/L)
    prop.err <- 0.15   #proportional error 
    # Initial estimates of the variances of the inter-individual variability
    eta.ka ~ 0.5   
    eta.cl ~ 0.1   
    eta.v  ~ 0.1   
  })
  model({
    # Where the model is specified
    # Individual parameters: typical value and inter-individual variability on the log scale
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv + eta.v)
    ka <- exp(lka + eta.ka)
    ## ODE example
    d/dt(depot)   = -ka * depot
    d/dt(central) =  ka * depot - (cl/v) * central
    ## Concentration is calculated
    cp = central/v
    ## And is assumed to follow proportional and additive error
    cp ~ prop(prop.err) + add(add.err)
  })
}


# Fit the base model with FOCEI; print = 0 switches off the printing of the estimation progress
# and NPDE and CWRES are added to the output table
# run100 <-
#  nlmixr2(One.comp.KA.ODE,    #the model definition
#          PKdata,                #the data set
#          est = "focei", control = foceiControl(print = 0),
#          table=list(npde=TRUE, cwres=TRUE))


# better load now
load("run100.Rdata")



# Create the ggPMX controller from the fit
# conts = continuous covariates, cats = categorical covariates,
# vpc = FALSE -> no VPC is created, is.draft = FALSE -> no DRAFT label in the plots
ctr2 <- pmx_nlmixr(run100, conts = c("WT","AGE"),cats=c("SEX","SPARSE"), vpc=FALSE,settings=pmx_settings(is.draft=FALSE))


# GOF plots of the base model:
# NPDE vs. population predictions
ctr2 %>% pmx_plot_npde_pred
# NPD vs. population predictions
ctr2 %>% pmx_plot_npd_pred
# CWRES vs. population predictions
ctr2 %>% pmx_plot_cwres_pred
# Individual random effects (etas) vs. the continuous covariates
ctr2 %>% pmx_plot_eta_conts





# ---- 2. Alternative absorption model: one transit compartment -----------------

One.comp.transit <- function() {
  ini({
    # Where initial estimates are specified
    lktr <- log(1.15)  #log k transit (1/h)
    lcl  <- log(0.135) #log Cl (L/h)
    lv   <- log(8)     #log V (L)
    prop.err <- 0.15   #proportional error 
    add.err <- 0.6     #additive error (mg/L)
    # Initial estimates of the variances of the inter-individual variability
    eta.ktr ~ 0.5   
    eta.cl ~ 0.1   
    eta.v ~ 0.1  
  })
  model({
    # Individual parameters: typical value and inter-individual variability on the log scale
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv + eta.v)
    ktr <- exp(lktr + eta.ktr)
    # rxode2-style differential equations are supported
    # depot -> transit1 -> central, the same rate constant ktr is used for both transfer steps
    d/dt(depot)   = -ktr * depot
    d/dt(central) =  ktr * transit1 - (cl/v) * central
    d/dt(transit1)   =  ktr * depot - ktr * transit1
    ## Concentration is calculated
    cp = central/v
    # And is assumed to follow proportional and additive error
    cp ~ prop(prop.err) + add(add.err)
  })
}


# Fit the transit compartment model with the same settings as the base model
# run101 <-
#  nlmixr2(One.comp.transit,    #the model definition
#          PKdata,                #the data set
#          est = "focei", control = foceiControl(print = 0),
#          table=list(npde=TRUE, cwres=TRUE))


# better load now
load("run100.Rdata")

# Create the ggPMX controller from the fit (same settings as for the base model)
ctr3 <- pmx_nlmixr(run101, conts = c("WT","AGE"),cats=c("SEX","SPARSE"), vpc=FALSE,settings=pmx_settings(is.draft=FALSE))


# GOF plots of the transit compartment model (same plots as for the base model)
ctr3 %>% pmx_plot_npde_pred
ctr3 %>% pmx_plot_npd_pred
ctr3 %>% pmx_plot_cwres_pred
ctr3 %>% pmx_plot_eta_conts


#################################################################################
##                                                                             ##
## Implement allometric covariates                                             ##
##                                                                             ##
#################################################################################

## One compartment transit model with allometric scaling on WT
## (WT is normalised to 70 in log(WT/70); the exponents are fixed)

run101_allo  <- run101 |>
  model( cl <- exp(lcl + eta.cl + ALLC * log(WT/70))) |> # add allometric scaling on cl
  model(     v  <- exp(lv + eta.v + ALLV * log(WT/70)))|> # add allometric scaling on volume
  ini(ALLC = fix(0.75))|>#allometric exponent cl (fixed to 0.75)
  ini(ALLV = fix(1.00)) #allometric exponent v (fixed to 1)

# Fit the allometric model (print the estimation progress every 5th iteration)
run107 <-
  nlmixr(run101_allo ,
         PKdata,
         est = "focei",
         foceiControl(print = 5))
run107

## do you get a significant drop in OFV by including allometric weight?
# load a previously saved result from file (the file has to exist in the working directory;
# it is not saved in this script)

save(run107, file="run107.Rdata")
load(file="run107.Rdata")
# difference in OFV between the allometric model and run101
run107$OBJF-run101$OBJF


#################################################################################
##                                                                             ##
##  Hands-on assignments: nlmixr2 model development                            ##
##                                                                             ##
##  Exercise: How about free fitting the exponents?                            ##
##  Run the allometric model (run101_allo) without fixing the exponents        ##
##  ALLC and ALLV. Do you get a better fit?                                    ##
##                                                                             ##
##  Hint: unfix() releases the fixed initial estimates, so that they are       ##
##  estimated                                                                  ##
##                                                                             ##
##  Compare the OFV of your model with the OFV of run107 (fixed exponents)     ##
##                                                                             ##
#################################################################################

# Your code here
