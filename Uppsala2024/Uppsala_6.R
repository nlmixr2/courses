## load the required libraries
library(xpose.nlmixr2)
library(nlmixr2)
library(lattice)
library(data.table)
library(xgxr)
library(patchwork)

## specify array of colours for curves
nlmixCOLS <- c("#28466A","#8DB6CD","#B40000")  

#################################################################################
##                                                                             ##
## nlmixr2 analysis Part 6                                                      ##
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
load(file = "run007F.Rdata")

## Extract EBEs from the model object:
EBEs <- as.data.table(run007FF)
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

#Only one of these is required but doing both will not harm
PKPDdata[, CMT := ifelse(AMT > 0, "depot", ifelse(DVID == 1, "central","effect"))]
PKPDdata[,DVID:=as.character(factor(DVID,levels=c(1,2),labels=c("central","effect")))]

## merge the data files with the EBEs
PDdata <- merge(PKPDdata, EBEs, by = "ID", all.x = TRUE)

# ## remove PK measurements from file so as to not confuse the vpcPlot later on
PDdata[, DEL := ifelse(DVID == "central" & AMT == 0, 1, 0)]
PDdata<-PDdata[DEL==0]



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



ggpd1<-ggplot(data = PDdata[MDV==0], aes(x = TIME, y = DV)) +
  geom_line(aes(group = ID), color = "grey50", size = 1, alpha = 0.3) +
  xgx_geom_ci(aes(x = TIME, color = NULL, group = NULL, shape = NULL), conf_level = 0.95) +
  labs(y = "TOF Response (PCA, %)", x="Time (hours)",color = "Group") +
  theme(legend.position = "none") 
print(ggpd1)


## run nlmixr2
run008F <-
  nlmixr2(KA1tr1IPP_PDtoemax1, 
         PDdata,  
         est = "focei", 
         foceiControl(print = 5))
run008F
save(run008F, file = "run008F.Rdata")


## generate the VPC

vpc8F<-vpcPlot(
  run008F,
  n = 500,
  show = list(obs_dv = TRUE),    #additional items to show, like the observations
  xlab = "Time (h)",
  ylab = "TOF Response (PCA, %)"
)
vpc8F


## generate smooth individual profiles

indivs <- augPred(run008F)

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
## nlmixr2 analysis Part 7                                                      ##
##                                                                             ##
## Simultaneous PKPD analysis                                                  ##
##                                                                             ##
#################################################################################



#################################################################################
##                                                                             ##
## Immediate effect using Emax model                                           ##
##                                                                             ##
#################################################################################

KA1tr1_PDimmemax1 <- function() {
  ini({
    ## PK
    tktr <- log(1.3)   # log ktr (/h)
    tcl  <- log(0.15) # log CL (L/h)
    tv   <- log(8)   # log Vc (L)

    eta.ktr ~ 1
    eta.cl ~ 0.3
    eta.v ~ 0.3
    eps.pkprop <- 0.1
    eps.pkadd <- 0.4

    ## PD
    tc50  <- log(1.5)    #log ec50 (mg/L)
    te0   <- log(100)  #log e0

    eta.c50  ~ .5
    eta.e0 ~ .1

    eps.pdadd <- 20

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


nlmixr2(KA1tr1_PDimmemax1) # Show initial estimates and model

run009F <-
  nlmixr2(KA1tr1_PDimmemax1,
         PKPDdata,
         est = "focei",
         foceiControl(print = 5))
run009F
save(run009F, file = "run009F.Rdata")
#load(file = "run009F.Rdata")

# ## generate the VPC of both endpoints (currently broken code)
# 
# vpc9F <-
#   vpcPlot(
#     run009F,
#     n = 500,
#     show = list(obs_dv = TRUE),
#     title = "Immediate effect simultaneous PKPD model",
#     xlab = "Time (h)",
#     ylab = "Warfarin (mg/L)\nTOF Response (PCA, %)"
#   )
# vpc9F



#Using direct Ron Keizer code to allow separate panels from central and effect:
library(vpc)

SimVPC9F<-vpcSim(
  run009F,                  #the nlmixr object
  n = 500)


setnames(SimVPC9F,c("id","time","sim"),c("ID","TIME","DV"))

vpc9F<-vpc_vpc(sim = SimVPC9F, obs = PKPDdata,#bins=bin_edges,
        show = list(obs_dv = TRUE,
                    obs_median=TRUE,
                    sim_median=TRUE,
                    sim_median_ci=TRUE,
                    obs_ci=TRUE,
                    pi=TRUE),
        xlab = "Time (h)",              #x-axis label
        ylab = "Warfarin (mg/L)\nTOF Response (PCA, %)",  #y-axis label
        stratify="CMT",
        title = "Immediate effect simultaneous PKPD model")
vpc9F


indivs1 <- augPred(run009F)

## default plot
plot(indivs1)


## nicer plot with lattice
## effect:
indivs1<-as.data.table(indivs1)
xyplot(
  values ~ time/24 | id,
  data = indivs1[Endpoint=="effect"],
  groups = ind,
  type = c("l", "l", "p"),
  col = nlmixCOLS[c(2, 1, 3)],
  cex = c(0.1, 0.1, 1),
  lwd = c(2, 2, 2),
  pch = c(1, 1, 19),
  xlab = "Time (days)\n",
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

## concentration:
xyplot(
  values ~ time/24 | id,
  data = indivs1[Endpoint=="central"],
  groups = ind,
  type = c("l", "l", "p"),
  col = nlmixCOLS[c(2, 1, 3)],
  cex = c(0.1, 0.1, 1),
  lwd = c(2, 2, 2),
  pch = c(1, 1, 19),
  xlab = "Time (days)\n",
  ylab = "Warfarin (mg/L)",
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





#If you have an older version of nlmixr2, then augPred will not work in this specific case
#This is general code that can also be used, 
#but will only provide individual predicted curves not population predictions
#
#Generate a data.table of EBEs
# EBEs009F <- data.table(run009F)
# EBEs009F <-
#   EBEs009F[!duplicated(ID), .(ID, ktr, cl, v, c50, e0)]
# #And a data.table of doses
# ## Read in the PKPD data file:
# dataF <- fread("warfarin_dat.csv")
# setnames(dataF, names(dataF), toupper(names(dataF)))
# Doses <- dataF[AMT > 0, .(ID, TIME, AMT, EVID = 1)]
# #Generate an eventTable by combining the doses with the same fine-meshed sampling points for all subjects
# evt <- et(Doses) %>% et(0, 150, by = 0.5)
# 
# #Define the RxODE model (emax is fixed to 1)
# modPKPD1 <- RxODE({
#   cp           =  central/v
#   d/dt(depot)  = -ktr * depot
#   d/dt(trans)  =  ktr * depot - ktr * trans
#   d/dt(central)=  ktr * trans - cl * cp
#   effect       =  e0 * (1 - cp/(c50 + cp))
# })
# 
# #Solve the system:
# res1 <- data.table(rxSolve(modPKPD1, EBEs009F, evt))
# #Merge with observed data points
# Rdata <- PKPDdata[MDV == 0, .(id = ID, time = TIME, DVID, DV)]
# xx1 <- merge(res1, Rdata, by = c("id", "time"), all.x = TRUE)
# xx1[,ID:=factor(id)]
# 
# #And plot:
# xyplot(
#   DV + cp ~ time | ID,
#   data = xx1[DVID == 'central' | is.na(DVID)],
#   type = c("b", "l"),
#   col = nlmixCOLS[c(3, 1)],
#   main = "Warfarin profiles for simultaneous PKPD effect-compartment model with Emax fixed to 1",
#   cex = c(1, 0.1),
#   lwd = 2,
#   pch = c(19, 1),
#   xlab = "Time (h)\n",
#   ylab = "Warfarin (mg/L)",
#   as.table = TRUE,
#   layout=c(8,4),
#   scales = list(alternating = 1)
# )
# 
# xyplot(
#   DV + effect ~ time | ID,
#   data = xx1[DVID == 'effect' | is.na(DVID)],
#   type = c("b", "l"),
#   col = nlmixCOLS[c(3, 2)],
#   main = "PCA profiles for simultaneous PKPD effect-compartment model with Emax fixed to 1",
#   cex = c(1, 0.1),
#   lwd = 2,
#   pch = c(19, 1),
#   xlab = "Time (h)\n",
#   ylab = "TOF Response (PCA, %)",
#   as.table = TRUE,
#   layout=c(8,4),
#   scales = list(alternating = 1)
# )




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

nlmixr2(KA1tr1_PDceemax) # Show initial estimates and model

## run nlmixr2

run010F <-
  nlmixr2(KA1tr1_PDceemax, 
         PKPDdata, 
         est = "focei", 
         foceiControl(print = 5))
run010F
save(run010F, file = "run010F.Rdata")
#load(file="run010F.Rdata")

## generate the VPC
#
# vpc10F <-
#   vpcPlot(
#     run010F,
#     n = 500,
#     show = list(obs_dv = TRUE),
#     title = "Effect compartment simultaneous PKPD model",
#     xlab = "Time (h)",
#     ylab = "Warfarin (mg/L)\nTOF Response (PCA, %)"
#   )
# vpc10F



SimVPC10<-vpcSim(
  run010F,                  #the nlmixr object
  n = 500)


setnames(SimVPC10,c("id","time","sim"),c("ID","TIME","DV"))

vpc10F <-vpc_vpc(sim = SimVPC10, obs = PKPDdata,#bins=bin_edges,
        show = list(obs_dv = TRUE,
                    obs_median=TRUE,
                    sim_median=TRUE,
                    sim_median_ci=TRUE,
                    obs_ci=TRUE,
                    pi=TRUE),
        xlab = "Time (h)",              #x-axis label
        ylab = "Warfarin (mg/L)\nTOF Response (PCA, %)",  #y-axis label
        stratify="CMT",
        title = "Effect compartment simultaneous PKPD model")

vpc10F


#Generate individual profiles using augPred:

indivs2 <- augPred(run010F)

## default plot
plot(indivs2)



## nicer plot with lattice
## effect:
indivs2 <- as.data.table(indivs2)

xyplot(
  values ~ time/24 | id,
  data = indivs2[Endpoint=="effect"],
  groups = ind,
  type = c("l", "l", "p"),
  col = nlmixCOLS[c(2, 1, 3)],
  cex = c(0.1, 0.1, 1),
  lwd = c(2, 2, 2),
  pch = c(1, 1, 19),
  xlab = "Time (days)\n",
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

## concentration:
xyplot(
  values ~ time/24 | id,
  data = indivs2[Endpoint=="central"],
  groups = ind,
  type = c("l", "l", "p"),
  col = nlmixCOLS[c(2, 1, 3)],
  cex = c(0.1, 0.1, 1),
  lwd = c(2, 2, 2),
  pch = c(1, 1, 19),
  xlab = "Time (days)\n",
  ylab = "Warfarin (mg/L)",
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





#If you have an older version of nlmixr2, then augPred will not work in this specific case
#This is general code that can also be used,
#but will only provide individual predicted curves not population predictions
#
#Generate a data.table of EBEs
# EBEs010F <- data.table(run010F)
# EBEs010F <-
#   EBEs010F[!duplicated(ID), .(ID, ktr, cl, v, c50, kout, e0)]
# #And a data.table of doses
# ## Read in the PKPD data file:
# dataF <- fread("warfarin_dat.csv")
# setnames(dataF, names(dataF), toupper(names(dataF)))
# Doses <- dataF[AMT > 0, .(ID, TIME, AMT, EVID = 1)]
# #Generate an eventTable by combining the doses with the same fine-meshed sampling points for all subjects
# evt <- et(Doses) %>% et(0, 150, by = 0.5)
# 
# #Define the RxODE model (emax is fixed to 1)
# modPKPD2 <- RxODE({
#   cp           =  central/v
#   d/dt(depot)  = -ktr * depot
#   d/dt(trans)  =  ktr * depot - ktr * trans
#   d/dt(central)=  ktr * trans - cl * cp
#   d/dt(ce)     =  kout * (cp - ce)
#   effect       =  e0 * (1 - ce/(c50 + ce))
# })
# 
# #Solve the system:
# res2 <- data.table(rxSolve(modPKPD2, EBEs010F, evt))
# #Merge with observed data points
# Rdata <- PKPDdata[MDV == 0, .(id = ID, time = TIME, DVID, DV)]
# xx2 <- merge(res2, Rdata, by = c("id", "time"), all.x = TRUE)
# xx2[,ID:=factor(id)]
# 
# #And plot:
# xyplot(
#   DV + cp ~ time | ID,
#   data = xx2[DVID == 'central' | is.na(DVID)],
#   type = c("b", "l"),
#   col = nlmixCOLS[c(3, 1)],
#   main = "Warfarin profiles for simultaneous PKPD effect-compartment model with Emax fixed to 1",
#   cex = c(1, 0.1),
#   lwd = 2,
#   pch = c(19, 1),
#   xlab = "Time (h)\n",
#   ylab = "Warfarin (mg/L)",
#   as.table = TRUE,
#   layout=c(8,4),
#   scales = list(alternating = 1)
# )
# 
# xyplot(
#   DV + effect ~ time | ID,
#   data = xx2[DVID == 'effect' | is.na(DVID)],
#   type = c("b", "l"),
#   col = nlmixCOLS[c(3, 2)],
#   main = "PCA profiles for simultaneous PKPD effect-compartment model with Emax fixed to 1",
#   cex = c(1, 0.1),
#   lwd = 2,
#   pch = c(19, 1),
#   xlab = "Time (h)\n",
#   ylab = "TOF Response (PCA, %)",
#   as.table = TRUE,
#   layout=c(8,4),
#   scales = list(alternating = 1)
# )




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

nlmixr2(KA1tr1IPP_PDtoemax) # Show initial estimates and model

## run nlmixr

run011F <-
  nlmixr(KA1tr1IPP_PDtoemax, 
         PKPDdata, 
         est = "focei", 
         foceiControl(print = 5))
run011F
save(run011F, file = "run011F.Rdata")
#load(file="run011F.Rdata")

## generate the VPC

# vpc11F <-
#   vpcPlot(
#     run011F,
#     n = 500,
#     show = list(obs_dv = TRUE),
#     title = "Turnover simultaneous PKPD model",
#     xlab = "Time (h)",
#     ylab = "Warfarin (mg/L)\nTOF Response (PCA, %)"
#   )
# vpc11F

SimVPC11<-vpcSim(
  run011F,                  #the nlmixr object
  n = 500)


setnames(SimVPC11,c("id","time","sim"),c("ID","TIME","DV"))

vpc11F <-vpc_vpc(sim = SimVPC11, obs = PKPDdata,#bins=bin_edges,
                 show = list(obs_dv = TRUE,
                             obs_median=TRUE,
                             sim_median=TRUE,
                             sim_median_ci=TRUE,
                             obs_ci=TRUE,
                             pi=TRUE),
                 xlab = "Time (h)",              #x-axis label
                 ylab = "Warfarin (mg/L)\nTOF Response (PCA, %)",  #y-axis label
                 stratify="CMT",
                 title = "Turnover simultaneous PKPD model")

vpc11F




## using patchwork syntax
vpc9F/vpc10F/vpc11F


indivs3 <- augPred(run011F)

## default plot
plot(indivs3)

## nicer plot with lattice
## effect:
indivs3<-as.data.table(indivs3)
xyplot(
  values ~ time/24 | id,
  data = indivs3[Endpoint=="effect"],
  groups = ind,
  type = c("l", "l", "p"),
  col = nlmixCOLS[c(2, 1, 3)],
  cex = c(0.1, 0.1, 1),
  lwd = c(2, 2, 2),
  pch = c(1, 1, 19),
  xlab = "Time (days)\n",
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

## concentration:
xyplot(
  values ~ time/24 | id,
  data = indivs3[Endpoint=="central"],
  groups = ind,
  type = c("l", "l", "p"),
  col = nlmixCOLS[c(2, 1, 3)],
  cex = c(0.1, 0.1, 1),
  lwd = c(2, 2, 2),
  pch = c(1, 1, 19),
  xlab = "Time (days)\n",
  ylab = "Warfarin (mg/L)",
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





