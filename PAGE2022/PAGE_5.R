## load the required libraries
library(xpose.nlmixr2)
library(nlmixr2)
library(lattice)
library(data.table)
library(patchwork)

## specify array of colours for curves
nlmixCOLS <- c("#28466A","#8DB6CD","#B40000")  

#################################################################################
##                                                                             ##
## nlmixr analysis Part 6                                                      ##
##                                                                             ##
## Sequential PKPD analysis                                                    ##
##                                                                             ##
#################################################################################



#################################################################################
##                                                                             ##
## Generate a data file that combines the PD measurements                      ##
## with the PK EBEs and PK doses                                               ##
##                                                                             ##
#################################################################################

## Use the results from the one compartment transit model to start the PD analysis
load(file = "fitOne.comp.transit.allo_F.Rdata")

## Extract EBEs from the model object:
EBEs <- as.data.table(fitOne.comp.transit.allo_F)
EBEs <-
  EBEs[!duplicated(ID), .(
    ID = as.numeric(as.character(ID)),
    IKTR = ktr,
    ICL = cl,
    IV = v
  )]

## Read in the PKPD data file:
PKPDdata <- fread("warfarin_dat.csv")
setnames(PKPDdata, names(PKPDdata), toupper(names(PKPDdata)))
## define MDV data items
PKPDdata[, MDV := ifelse(is.na(DV), 1, 0)]
PKPDdata[, MDV := ifelse(AMT > 0, 1, MDV)]

## merge the data files with the EBEs
PDdata <- merge(PKPDdata, EBEs, by = "ID", all.x = TRUE)

## mark the PK measurements from the PKPD file to not be analysed by setting MDV to 1
PDdata[, MDV := ifelse(DVID == 1 & AMT == 0, 1, MDV)]
PDdata[,DVID:=NULL]

#################################################################################
##                                                                             ##
## Define a turnover model and use EBEs for PK                                 ##
##                                                                             ##
#################################################################################


KA1tr1IPP_PDtoemax1 <- function() {
  ini({
    tc50  <- log(1)    #log ec50 (mg/L)
    tkout <- log(0.05) #log tkout (/h)
    te0   <- log(100)  #log e0
    
    eta.c50  ~ .5
    eta.kout ~ .1
    eta.e0 ~ .1
    
    eps.pdadd <- 100
  })
  model({
    c50 =  exp(tc50 + eta.c50)
    kout = exp(tkout + eta.kout)
    e0 = exp(te0 + eta.e0)
    
    # PK parameters from input dataset
    ktr = IKTR
    cl = ICL
    v = IV
    
    cp = central/v
    
    PD = 1 - cp/(c50 + cp)
    
    effect(0) = e0
    kin = e0 * kout
    
    d/dt(depot)   = -ktr * depot
    d/dt(central) =  ktr * trans - cl/v * central
    d/dt(trans)   =  ktr * depot - ktr * trans
    
    d/dt(effect)  =  kin * PD - kout * effect
    
    effect ~ add(eps.pdadd)
  })
}

## run nlmixr

fitKA1tr1IPP_PDtoemax1_F <-
  nlmixr(KA1tr1IPP_PDtoemax1, 
         PDdata,  
         est = "focei", 
         foceiControl(print = 5))
fitKA1tr1IPP_PDtoemax1_F
save(fitKA1tr1IPP_PDtoemax1_F, file = "fitKA1tr1IPP_PDtoemax1_F.Rdata")

## generate the VPC

vpcPlot(
  fitKA1tr1IPP_PDtoemax1_F,
  n = 500,
  show = list(obs_dv = TRUE,
              obs_median=TRUE,
              sim_median=TRUE,
              sim_median_ci=TRUE,
              obs_ci=TRUE,pi=TRUE),    #additional items to show, like the observations
  xlab = "Time (h)",
  ylab = "TOF Response (PCA, %)"
)

## generate smooth individual profiles

indivs <- augPred(fitKA1tr1IPP_PDtoemax1_F)

## default plot
plot(indivs)

## nicer plot with lattice
xyplot(
  values ~ time | id,
  data = indivs,
  groups = ind,
  type = c("l", "l", "p"),
  col = nlmixCOLS[c(2, 1, 3)],
  cex = c(0.1, 0.1, 1),
  lwd = c(2, 2, 2),
  pch = c(1, 1, 19),
  xlab = "Time (h)\n",
  ylab = "TOF Response (PCA, %)",
  layout=c(8,4),
  as.table = TRUE,
  scales = list(alternating = 1),
  auto.key = list(
    adj = 1,
    col = nlmixCOLS[c(2, 1, 3)],
    columns = 3,
    space = "bottom",
    rectangles = FALSE,
    points = FALSE
  )
)


#################################################################################
##                                                                             ##
## nlmixr analysis Part 7                                                      ##
##                                                                             ##
## Simultaneous PKPD analysis                                                  ##
##                                                                             ##
#################################################################################


PKPDdata[,dvid:=as.character(factor(DVID,levels=c(1,2),labels=c("central","effect")))]
setnames(PKPDdata,"DVID","DVIDold")

#################################################################################
##                                                                             ##
## Immediate effect using Emax model                                           ##
##                                                                             ##
#################################################################################

KA1tr1_PDimmemax1 <- function() {
  ini({
    ## PK
    tktr <- log(1)   # log ktr (/h)
    tcl  <- log(0.1) # log CL (L/h)
    tv   <- log(8)   # log Vc (L)
    
    eta.ktr ~ 1
    eta.cl ~ 0.1
    eta.v ~ 0.1
    eps.pkprop <- 0.1
    eps.pkadd <- 0.4
    
    ## PD
    tc50  <- log(1)    #log ec50 (mg/L)
    te0   <- log(100)  #log e0
    
    eta.c50  ~ .5
    eta.e0 ~ .1
    
    eps.pdadd <- 100
    
  })
  model({
    ktr <- exp(tktr + eta.ktr)
    cl  <- exp(tcl + eta.cl)
    v   <- exp(tv + eta.v)
    
    c50  = exp(tc50 + eta.c50)
    e0   = exp(te0 + eta.e0)
    
    cp           =  central/v
    d/dt(depot)  = -ktr * depot
    d/dt(central)=  ktr * trans - cl * cp
    d/dt(trans)  =  ktr * depot - ktr * trans
    effect       =  e0 * (1 - cp/(c50 + cp))
    
    cp ~ prop(eps.pkprop) + add(eps.pkadd) | central
    effect ~ add(eps.pdadd) | effect
  })
}

nlmixr(KA1tr1_PDimmemax1) # Show initial estimates and model

fitKA1tr1_PDimmemax1_F <-
  nlmixr(KA1tr1_PDimmemax1,
         PKPDdata, 
         est = "focei", 
         foceiControl(print = 5))
fitKA1tr1_PDimmemax1_F
save(fitKA1tr1_PDimmemax1_F, file = "fitKA1tr1_PDimmemax1_F.Rdata")
#load(file = "fitKA1tr1_PDimmemax1_F.Rdata")

## generate the VPC of both endpoints

vpcKA1tr1_PDimmemax1 <-
  vpcPlot(
    fitKA1tr1_PDimmemax1_F,
    n = 500,
    show = list(obs_dv = TRUE),
    title = "Immediate effect simultaneous PKPD model",
    xlab = "Time (h)",
    ylab = "Warfarin (mg/L)\nTOF Response (PCA, %)"
  )
vpcKA1tr1_PDimmemax1



#################################################################################
##                                                                             ##
## Effect compartment model using Emax                                         ##
##                                                                             ##
#################################################################################



KA1tr1_PDceemax <- function() {
  ini({
    ## PK
    tktr <- log(1)   # log ktr (/h)
    tcl  <- log(0.1) # log CL (L/h)
    tv   <- log(8)   # log Vc (L)
    
    eta.ktr ~ 1
    eta.cl ~ 0.1
    eta.v ~ 0.1
    eps.pkprop <- 0.1
    eps.pkadd <- 0.4
    
    ## PD
    tc50  <- log(1)    #log ec50 (mg/L)
    tkout <- log(0.05) #log tkout (/h)
    te0   <- log(100)  #log e0
    
    eta.c50  ~ .5
    eta.kout ~ .1
    eta.e0 ~ .1
    
    eps.pdadd <- 100
  })
  model({
    ktr <- exp(tktr + eta.ktr)
    cl  <- exp(tcl + eta.cl)
    v   <- exp(tv + eta.v)
    
    c50  = exp(tc50 + eta.c50)
    kout = exp(tkout + eta.kout)
    e0   = exp(te0 + eta.e0)
    emax = 1
    
    cp           =  central/v
    d/dt(depot)  = -ktr * depot
    d/dt(central) =  ktr * trans - cl * cp
    d/dt(trans)  =  ktr * depot - ktr * trans
    d/dt(ce)     =  kout * (cp - ce)
    
    effect       =  e0 * (1 - emax * ce/(c50 + ce))
    
    cp ~ prop(eps.pkprop) + add(eps.pkadd) | central
    effect ~ add(eps.pdadd) | effect
  })
}

nlmixr(KA1tr1_PDceemax) # Show initial estimates and model

## run nlmixr

fitKA1tr1_PDceemax_F <-
  nlmixr(KA1tr1_PDceemax, 
         PKPDdata, 
         est = "focei", 
         foceiControl(print = 5))
fitKA1tr1_PDceemax_F
save(fitKA1tr1_PDceemax_F, file = "fitKA1tr1_PDceemax_F.Rdata")
#load(file="fitKA1tr1_PDceemax_F.Rdata")

## generate the VPC

vpcKA1tr1_PDceemax <-
  vpcPlot(
    fitKA1tr1_PDceemax_F,
    n = 500,
    show = list(obs_dv = TRUE),
    title = "Effect compartment simultaneous PKPD model",
    xlab = "Time (h)",
    ylab = "Warfarin (mg/L)\nTOF Response (PCA, %)"
  )
vpcKA1tr1_PDceemax



#Generate a data.table of EBEs
EBEs_KA1tr1_PDceemax <- data.table(fitKA1tr1_PDceemax_F)
EBEs_KA1tr1_PDceemax <-
  EBEs_KA1tr1_PDceemax[!duplicated(ID), .(ID, ktr, cl, v, c50, kout, e0)]
#And a data.table of doses
Doses <- dataF[AMT > 0, .(ID, TIME, AMT, EVID = 1)]
#Generate an eventTable by combining the doses with the same fine-meshed sampling points for all subjects
evt <- et(Doses) %>% et(0, 150, by = 0.5)

#Define the RxODE model
modPKPD1 <- RxODE({
  cp           =  central/v
  d/dt(depot)  = -ktr * depot
  d/dt(trans)  =  ktr * depot - ktr * trans
  d/dt(central)=  ktr * trans - cl * cp
  d/dt(ce)     =  kout * (cp - ce)
  effect       =  e0 * (1 - ce/(c50 + ce))
})

#Solve the system:
res1 <- data.table(rxSolve(modPKPD1, EBEs_KA1tr1_PDceemax, evt))
#Merge with observed data points
Rdata <- PKPDdata[MDV == 0, .(id = ID, time = TIME, dvid, DV)]
xx1 <- merge(res1, Rdata, by = c("id", "time"), all.x = TRUE)
xx1[,ID:=factor(id)]

#And plot:
xyplot(
  DV + cp ~ time | ID,
  data = xx1[dvid == 'central' | is.na(dvid)],
  type = c("b", "l"),
  col = nlmixCOLS[c(3, 1)],
  main = "Warfarin profiles for simultaneous PKPD effect-compartment model with Emax fixed to 1",
  cex = c(1, 0.1),
  lwd = 2,
  pch = c(19, 1),
  xlab = "Time (h)\n",
  ylab = "Warfarin (mg/L)",
  as.table = TRUE,
  layout=c(8,4),
  scales = list(alternating = 1)
)

xyplot(
  DV + effect ~ time | ID,
  data = xx1[dvid == 'effect' | is.na(dvid)],
  type = c("b", "l"),
  col = nlmixCOLS[c(3, 2)],
  main = "PCA profiles for simultaneous PKPD effect-compartment model with Emax fixed to 1",
  cex = c(1, 0.1),
  lwd = 2,
  pch = c(19, 1),
  xlab = "Time (h)\n",
  ylab = "TOF Response (PCA, %)",
  as.table = TRUE,
  layout=c(8,4),
  scales = list(alternating = 1)
)




#################################################################################
##                                                                             ##
## Turnover model using Emax                                                   ##
##                                                                             ##
#################################################################################


KA1tr1IPP_PDtoemax <- function() {
  ini({
    ## PK
    tktr <- log(1)   # log ktr (/h)
    tcl  <- log(0.1) # log CL (L/h)
    tv   <- log(8)   # log Vc (L)
    
    eta.ktr ~ 1
    eta.cl ~ 0.1
    eta.v ~ 0.1
    eps.pkprop <- 0.1
    eps.pkadd <- 0.4
    
    ## PD
    tc50  <- log(1)    #log ec50 (mg/L)
    tkout <- log(0.05) #log tkout (/h)
    te0   <- log(100)  #log e0
    
    eta.c50  ~ .5
    eta.kout ~ .1
    eta.e0 ~ .1
    
    eps.pdadd <- 100
    
  })
  model({
    ktr <- exp(tktr + eta.ktr)
    cl  <- exp(tcl + eta.cl)
    v   <- exp(tv + eta.v)
    
    c50  = exp(tc50 + eta.c50)
    kout = exp(tkout + eta.kout)
    e0   = exp(te0 + eta.e0)
    emax = 1
    
    cp           =  central/v
    d/dt(depot)  = -ktr * depot
    d/dt(central)=  ktr * trans - cl * cp
    d/dt(trans)  =  ktr * depot - ktr * trans
    effect(0)    =  e0
    kin          =  e0 * kout
    PD           =  1 - emax * cp/(c50 + cp)
    d/dt(effect) =  kin * PD - kout * effect
    
    cp ~ prop(eps.pkprop) + add(eps.pkadd) | central
    effect ~ add(eps.pdadd) | effect
  })
}

nlmixr(KA1tr1IPP_PDtoemax) # Show initial estimates and model

## run nlmixr

fitKA1tr1IPP_PDtoemax_F <-
  nlmixr(KA1tr1IPP_PDtoemax, 
         PKPDdata, 
         est = "focei", 
         foceiControl(print = 5))
fitKA1tr1IPP_PDtoemax_F
save(fitKA1tr1IPP_PDtoemax_F, file = "fitKA1tr1IPP_PDtoemax_F.Rdata")
#load(file="fitKA1tr1IPP_PDtoemax_F.Rdata")

## generate the VPC

vpcKA1tr1IPP_PDtoemax <-
  vpcPlot(
    fitKA1tr1IPP_PDtoemax_F,
    n = 500,
    show = list(obs_dv = TRUE),
    title = "Turnover simultaneous PKPD model",
    xlab = "Time (h)",
    ylab = "Warfarin (mg/L)\nTOF Response (PCA, %)"
  )
vpcKA1tr1IPP_PDtoemax



## using patchwork syntax

vpcKA1tr1_PDimmemax1/vpcKA1tr1_PDceemax/vpcKA1tr1IPP_PDtoemax






#################################################################################
##                                                                             ##
## RxODE code Part 8                                                           ##
##                                                                             ##
## Bonus code to generate individual profiles by simulating both endpoints     ##
##                                                                             ##
#################################################################################

## Read in the PKPD data file:
dataF <- fread("warfarin_dat.csv")
setnames(dataF, names(dataF), toupper(names(dataF)))
## Generate a data.table of doses
Doses <- dataF[AMT > 0, .(ID, TIME, AMT, EVID = 1)]
## Generate an eventTable by combining the doses with the same fine-meshed 
## sampling points for all subjects using magrittr piping (%>%)
evt <- et(Doses) %>% et(0, 150, by = 0.5)

## Generate a data.table of EBEs
EBEs_KA1tr1IPP_PDtoemax <- data.table(fitKA1tr1IPP_PDtoemax_F)
EBEs_KA1tr1IPP_PDtoemax <-
  EBEs_KA1tr1IPP_PDtoemax[!duplicated(ID), .(ID, ktr, cl, v, c50, kout, e0)]

## Define the RxODE model
modPKPD2 <- RxODE({
  cp           =  central/v
  d/dt(depot)  = -ktr * depot
  d/dt(central)=  ktr * trans - cl * cp
  d/dt(trans)  =  ktr * depot - ktr * trans
  effect(0)    =  e0
  kin          =  e0 * kout
  PD           =  1 - cp/(c50 + cp)
  d/dt(effect) =  kin * PD - kout * effect
})

## Solve the system:
res <- data.table(rxSolve(modPKPD2, EBEs_KA1tr1IPP_PDtoemax, evt))
## Merge with observed data points
Rdata <- PKPDdata[MDV == 0, .(id = ID, time = TIME, dvid, DV)]
xx <- merge(res, Rdata, by = c("id", "time"), all.x = TRUE)
xx[,ID:=factor(id)]
## And plot:
xyplot(
  DV + cp ~ time/24 |ID,
  data = xx[dvid == 'central' | is.na(dvid)],
  type = c("b", "l"),
  col = nlmixCOLS[c(3, 1)],
  main = "Warfarin profiles for simultaneous PKPD turnover model",
  cex = c(1, 0.1),
  lwd = 2,
  pch = c(19, 1),
  xlab = "Time (days)\n",
  ylab = "Warfarin (mg/L)",
  as.table = TRUE,
  layout=c(8,4),
  scales = list(alternating = 1)
)

xyplot(
  DV + effect ~ time/24 | ID,
  data = xx[dvid == 'effect' | is.na(dvid)],
  type = c("b", "l"),
  col = nlmixCOLS[c(3, 2)],
  main = "PCA profiles for simultaneous PKPD turnover model",
  cex = c(1, 0.1),
  lwd = 2,
  pch = c(19, 1),
  xlab = "Time (days)\n",
  ylab = "TOF Response (PCA, %)",
  as.table = TRUE,
  layout=c(8,4),
  scales = list(alternating = 1)
)

