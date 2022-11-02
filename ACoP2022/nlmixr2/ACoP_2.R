#################################################################################
##                                                                             ##
##  Solutions to hands-on assignments: nlmixr model development                ##
##  Examine the GOF plots and implement models with                            ##
##  one or five transit compartments                                           ##
##                                                                             ##
##  Compare vpcs of alternatives and compare OFVs:                             ##
##                                                                             ##
##  fitPK001$OBJF-fitPK002$OBJF                                                ##
##                                                                             ##
#################################################################################


## load the required libraries
library(nlmixr2)
library(ggPMX)
library(patchwork)
library(ggplot2)
library(xpose.nlmixr2)
library(tidyverse)

if (!exists("PKdata", globalenv())) {
  source("ACoP_1.R")
}


#################################################################################
##                                                                             ##
##   Update the model with a single transit compartment                        ##
##                                                                             ##
#################################################################################

## One compartment transit model
One.comp.transit <- function() {
  ini({
    # Where initial conditions/variables are specified
    lktr <- log(1.15)  #log k transit (/h)
    lcl  <- log(0.135) #log Cl (L/hr)
    lv   <- log(8)     #log V (L)
    prop.err <- 0.15   #proportional error (SD/mean)
    add.err <- 0.6     #additive error (mg/L)
    eta.ktr ~ 0.5   
    eta.cl ~ 0.1   
    eta.v ~ 0.1  
  })
  model({
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv + eta.v)
    ktr <- exp(lktr + eta.ktr)
    # RxODE-style differential equations are supported
    d/dt(depot)   = -ktr * depot
    d/dt(central) =  ktr * trans - (cl/v) * central
    d/dt(trans)   =  ktr * depot - ktr * trans
    ## Concentration is calculated
    cp = central/v
    # And is assumed to follow proportional and additive error
    cp ~ prop(prop.err) + add(add.err)
  })
}



fitOne.comp.transit_F <-
  nlmixr(One.comp.transit,
         PKdata,
         est = "focei")


#load(file="fitOne.comp.transit_F.Rdata")

p3<-vpcPlot( 
  fitOne.comp.transit_F,                  #the nlmixr object
  bins=bin_edges,
  n = 500,                       #number of replicates simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE,
              obs_median=TRUE,
              sim_median=TRUE,
              sim_median_ci=TRUE,
              obs_ci=TRUE,
              pi=TRUE),    #additional items to show, like the observations
  xlab = "Time (h)",             #x-axis label
  ylab = "Concentration (mg/L)", #y-axis label
  title= "One transit compartment"
)

print(p3)


#using patchwork syntax:
print(p1/p3)


## restrict x-axis to first 12 hours:
p3m<-p3+xlim(0,12) 
p1m<-p1+xlim(0,12)

#Comparing the VPCs for the first 12 hours to see if absorption is described better using a transit compartment
print(p1m/p3m)



xpdb1 <- xpose_data_nlmixr(fitOne.comp.KA.ODE_F)
IWRES1<-absval_res_vs_idv(xpdb1,        #the xpose object
                            res = "IWRES",  #examine absolute values (absval) of individual weighted residuals
                            idv = "TIME",   #as a function of time
                            caption = NULL,#if not NULL provides the directory where this was run
                            title="No transit compartments",
                            subtitle=NULL) 


print(IWRES1)

xpdb3 <- xpose_data_nlmixr(fitOne.comp.transit_F)

#Absolute values of individual weighted residual vs time
IWRES3<-absval_res_vs_idv(xpdb3,        #the xpose object
                          res = "IWRES",  #examine absolute values (absval) of individual weighted residuals
                          idv = "TIME",   #as a function of time
                          caption = NULL,#if not NULL provides the directory where this was run
                          title="One transit compartment",
                          subtitle=NULL) 


print(IWRES3)



print(IWRES1/IWRES3)

#################################################################################
##                                                                             ##
##   Update the model with five transit compartments                            ##
##                                                                             ##
#################################################################################

## 5 transit compartments
KA1tr5ode <- function() {
  ini({
    # Where initial conditions/variables are specified
    lktr  <- log(1.15) #log transit rate constant (/h)
    lcl  <- log(0.135) #log Cl (L/h)
    lv   <- log(8)     #log V (L)
    prop.err <- 0.15   #proportional error (SD/mean)
    add.err  <- 0.6    #additive error (mg)
    eta.ktr ~ 0.5                   
    eta.cl ~ 0.1    
    eta.v  ~ 0.1    
  })
  model({
    # Where the model is specified
    ktr <- exp(lktr + eta.ktr)
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv + eta.v)
    ## ODE example
    cc=central/v
    d/dt(depot)= - ktr*depot
    d/dt(central) = ktr*transit5 - cl*cc
    d/dt(transit1)= ktr*(depot - transit1) 
    d/dt(transit2)= ktr*(transit1 - transit2) 
    d/dt(transit3)= ktr*(transit2 - transit3) 
    d/dt(transit4)= ktr*(transit3 - transit4) 
    d/dt(transit5)= ktr*(transit4 - transit5) 
    
    ## where residual error is assumed to follow proportional and additive error
    cc ~ prop(prop.err) + add(add.err)
  })
}

nlmixr(KA1tr5ode)

#################################################################################
##                                                                             ##
##   Run using FOCEi                                                           ##
##                                                                             ##
#################################################################################

fitKA1tr5ode_F <-
  nlmixr(KA1tr5ode,
         PKdata,
         est = "focei")

fitKA1tr5ode_F

save(fitKA1tr5ode_F,file="fitKA1tr5ode_F.Rdata")
#load(file="fitKA1tr5ode_F.Rdata")

p5<-vpcPlot( 
  fitKA1tr5ode_F,                 #the nlmixr object
  bins=bin_edges,
  n = 500,                        #number of replicates simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE,
              obs_median=TRUE,
              sim_median=TRUE,
              sim_median_ci=TRUE,
              obs_ci=TRUE,
              pi=TRUE),    #additional items to show, like the observations
  xlab = "Time (h)",             #x-axis label
  ylab = "Concentration (mg/L)", #y-axis label
  title= "Five transit compartments"
)
print(p5)

p5m<-p5+xlim(0,12) 
p3m<-p3+xlim(0,12) 
p1m<-p1+xlim(0,12)

print(p1m/p3m/p5m)




xpdb5 <- xpose_data_nlmixr(fitKA1tr5ode_F)

IWRES5<-absval_res_vs_idv(xpdb5,        #the xpose object
                          res = "IWRES",  #examine absolute values (absval) of individual weighted residuals
                          idv = "TIME",   #as a function of time
                          caption = NULL,#if not NULL provides the directory where this was run
                          title="Five transit compartments",
                          subtitle=NULL) 


print(IWRES5)



print(IWRES1/IWRES3/IWRES5)
prnt(IWRES1+IWRES3+IWRES5)



ctr1 <- pmx_nlmixr(fitOne.comp.KA.ODE_F, conts = c("WT","AGE"),cats=c("SEX","SPARSE"), vpc=FALSE,settings=pmx_settings(is.draft=FALSE))
ctr3 <- pmx_nlmixr(fitOne.comp.transit_F, conts = c("WT","AGE"),cats=c("SEX","SPARSE"), vpc=FALSE,settings=pmx_settings(is.draft=FALSE))
ctr5 <- pmx_nlmixr(fitKA1tr5ode_F, conts = c("WT","AGE"),cats=c("SEX","SPARSE"), vpc=FALSE,settings=pmx_settings(is.draft=FALSE))

p1cwres<-pmx_plot_cwres_time(ctr1,labels = list(x = "Time (hours)", y = "CWRESI",
                                                title="No transit compartments FOCEi",subtitle="CWRESI vs Time"))
p1npd<-pmx_plot_npd_time(ctr1,labels = list(x = "Time (hours)", y = "Normalized PD",
                                            title="No transit compartments FOCEi",subtitle="NPD vs Time"))

p3cwres<-pmx_plot_cwres_time(ctr3,labels = list(x = "Time (hours)", y = "CWRESI",
                                                title="One transit compartment FOCEi",subtitle="CWRESI vs Time"))
p3npd<-pmx_plot_npd_time(ctr3,labels = list(x = "Time (hours)", y = "Normalized PD",
                                            title="One transit compartment FOCEi",subtitle="NPD vs Time"))


p5cwres<-pmx_plot_cwres_time(ctr5,labels = list(x = "Time (hours)", y = "CWRESI",
                                                title="Five transit compartments FOCEi",subtitle="CWRESI vs Time"))
p5npd<-pmx_plot_npd_time(ctr5,labels = list(x = "Time (hours)", y = "Normalized PD",
                                            title="Five transit compartments FOCEi",subtitle="NPD vs Time"))




print(p1cwres+p3cwres+p5cwres)


print(p1npd+p3npd+p5npd)


fitOne.comp.KA.ODE_F$OBJF

fitOne.comp.transit_F$OBJF

fitOne.comp.transit_F$OBJF-fitOne.comp.KA.ODE_F$OBJF

fitKA1tr5ode_F$OBJF


fitKA1tr5ode_F$OBJF - fitOne.comp.KA.ODE_F$OBJF

fitKA1tr5ode_F$OBJF - fitOne.comp.transit_F$OBJF


