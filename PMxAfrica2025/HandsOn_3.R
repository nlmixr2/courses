#################################################################################
##                                                                             ##
##  Solutions to hands-on assignments: nlmixr2 model development               ##
##  Examine the VPCs and implement models with                                 ##
##  one or five transit compartments                                           ##
##                                                                             ##
#################################################################################

## load the required libraries
library(nlmixr2)
library(data.table)
library(ggPMX)
library(patchwork)
library(ggplot2)

## read in the Warfarin PK-only data set using data.table syntax (fast and efficient!)
PKdata <- fread("warfarin_PKS.csv")


#################################################################################
##                                                                             ##
##   Read in the previous model without transit compartment                    ##
##                                                                             ##
#################################################################################

One.comp.KA.ODE <- function() {
  ini({
    # Where initial conditions/variables are specified
    lka  <- log(1.15)  #log ka (1/h)
    lcl  <- log(0.135) #log Cl (L/h)
    lv   <- log(8)     #log V (L)
    prop.err <- 0.15   #proportional error (SD/mean)
    add.err  <- 0.6    #additive error (mg/L)
    eta.ka ~ 0.5   
    eta.cl ~ 0.1   
    eta.v  ~ 0.1   
  })
  model({
    # Where the model is specified
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
load(file = "run002F.Rdata")


## because the data set uses nominal time points, it is nice to have the bins surround these time points 
## so that each time point falls in a bin

bin_mids <- sort(unique(PKdata$TIME))
bin_edges <- bin_mids - c(0, diff(bin_mids) / 2) # puts the edges in the middle of the nominal time points

p2<-vpcPlot( 
  run002F,                       #the nlmixr2 object
  bins=bin_edges,                #the edges of the bins
  n = 500,                       #number of replicates simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE,
              obs_median=TRUE,
              sim_median=TRUE,
              sim_median_ci=TRUE,
              obs_ci=TRUE,
              pi=TRUE),
  xlab = "Time (h)",             #x-axis label
  ylab = "Concentration (mg/L)", #y-axis label
  title= "No transit compartments"
)
p2





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



run003F <-
  nlmixr(One.comp.transit,
         PKdata,
         est = "focei")
save(run003F,file="run003F.Rdata")
run003F
#load(file="run003F.Rdata")

p3<-vpcPlot( 
  run003F,                       #the nlmixr2 object
  bins=bin_edges,                #the edges of the bins
  n = 500,                       #number of replicates simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE,
              obs_median=TRUE,
              sim_median=TRUE,
              sim_median_ci=TRUE,
              obs_ci=TRUE,
              pi=TRUE),
  xlab = "Time (h)",             #x-axis label
  ylab = "Concentration (mg/L)", #y-axis label
  title= "One transit compartment"
)
p3




#using patchwork syntax:
p2/p3


## restrict x-axis to first 12 hours:
p3m<-p3+xlim(0,12) 
p2m<-p2+xlim(0,12)

#Comparing the VPCs for the first 12 hours to see if absorption is described better using a transit compartment
p2m/p3m





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
nlmixr2(KA1tr5ode)

#################################################################################
##                                                                             ##
##   Run using FOCEi                                                           ##
##                                                                             ##
#################################################################################

run004F <-
  nlmixr2(KA1tr5ode,
         PKdata,
         est = "focei")
run004F
save(run004F,file="run004F.Rdata")
#load(file="run004F.Rdata")

p5<-vpcPlot( 
  run004F,                       #the nlmixr2 object
  bins=bin_edges,                #the edges of the bins
  n = 500,                       #number of replicates simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE,
              obs_median=TRUE,
              sim_median=TRUE,
              sim_median_ci=TRUE,
              obs_ci=TRUE,
              pi=TRUE),
  xlab = "Time (h)",             #x-axis label
  ylab = "Concentration (mg/L)", #y-axis label
  title= "Five transit compartments"
)
p5

p5m<-p5+xlim(0,12) 
p3m<-p3+xlim(0,12) 
p2m<-p2+xlim(0,12)

#Stacked
p2m/p3m/p5m

#Side by side
p2m+p3m+p5m




run002F$OBJF

run003F$OBJF

run003F$OBJF-run002F$OBJF

run004F$OBJF

run004F$OBJF - run002F$OBJF

run004F$OBJF - run003F$OBJF


