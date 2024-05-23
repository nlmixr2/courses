## load the required libraries
library(xpose.nlmixr2)
library(nlmixr2)
library(lattice)
library(data.table)

## specify array of colours for curves
nlmixCOLS <- c("#28466A","#8DB6CD","#B40000")  

#################################################################################
##                                                                             ##
## nlmixr2 analysis Part 3                                                      ##
##                                                                             ##
## Generating Bayesian feedback estimates                                      ##
##                                                                             ##
#################################################################################

## read in the Warfarin PK-only data set using data.table syntax (fast and efficient!)
PKdata <- fread("warfarin_PKS.csv")


## Generating Bayesian feedback estimates
## nlmixr can generate empirical Bayes estimates for Bayesian feedback:
##  individual EBEs for a new data set using existing population parameters
## Useful in a therapeutic drug monitoring setting 
## Or for generating exposure estimates with a particularly nasty model that you 
##  do not want to refit on new data :-)

KA1tr5posthoc <- function() {
  ini({
    # Where initial conditions/variables are specified
    lktr <-  1.18994619   #log ktr (/h)
    lcl  <- -2.01737477   #log Cl (L/h)
    lv   <-  2.06631620   #log V (L)
    prop.err <- 0.07883633#proportional error (SD/mean)
    add.err <- 0.37249666 #additive error (mg/L)
    eta.ktr ~ 0.2532964   #IIV ktr
    eta.cl ~ 0.08073339   #IIV Cl
    eta.v ~ 0.04490733   #IIV V
  })
  model({
    # Where the model is specified
    # The model uses the ini-defined variable names
    cl  <- exp(lcl + eta.cl)
    v   <- exp(lv + eta.v)
    ktr <- exp(lktr + eta.ktr)
    # RxODE-style differential equations are supported
    d/dt(trns1) = -ktr * trns1
    d/dt(trns2) =  ktr * trns1 - ktr * trns2
    d/dt(trns3) =  ktr * trns2 - ktr * trns3
    d/dt(trns4) =  ktr * trns3 - ktr * trns4
    d/dt(trns5) =  ktr * trns4 - ktr * trns5
    d/dt(central) =  ktr * trns5 - (cl/v) * central
    ## Concentration is calculated
    cp = central/v
    # And is assumed to follow proportional and additive error
    cp ~ prop(prop.err) + add(add.err)
  })
}

## specify posthoc as estimation method to only generate EBEs and not do any estimation

run005ph <- nlmixr2(KA1tr5posthoc, PKdata, est = "posthoc")

##Generate individual smooth profiles using augPred
indivpk2 <- augPred(run005ph)

##and plot using lattice

xyplot(
  values ~ time | id,
  data = indivpk2,
  groups = ind,
  type = c("l", "l", "p"),
  col = nlmixCOLS[c(2, 1, 3)],
  cex = c(0.1, 0.1, 1),
  lwd = c(2, 2, 2),
  pch = c(1, 1, 19),
  xlab = "Time (h)\n",
  ylab = "Warfarin (mg/L)",
  layout=c(8,4),
  as.table = TRUE,
  scales = list(alternating = 1),
  main = "Five transit compartment model estimated using Bayesian feedback",
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
## nlmixr2 analysis Part 4                                                      ##
##                                                                             ##
## Implement mu-referenced covariates on log-scale                             ##
##                                                                             ##
#################################################################################


Covs <- PKdata[!duplicated(ID)]
table(Covs$SEX)
# 0  1
# 5 27

## One transit compartment model with Sex on V
KAtr1_sexV <- function() {
  ini({
    # Where initial conditions/variables are specified
    lktr <- log(1.15)  #log k transit (/h)
    lcl  <- log(0.135) #log CL (L/h)
    lv   <- log(8)     #log V (L)
    Sex_V <- 0.1       #log Sex on v
    prop.err <- 0.15   #proportional error (SD/mean)
    add.err <- 0.6     #additive error (mg/L)
    eta.ktr ~ 0.5   #IIV ktr
    eta.cl ~ 0.1   #IIV Cl
    eta.v ~ 0.1   #IIV V
  })
  model({
    #Sex on volume
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv + eta.v + Sex_V * SEX)
    ktr <- exp(lktr + eta.ktr)
    # RxODE-style differential equations are supported
    d/dt(depot) = -ktr * depot
    d/dt(central) =  ktr * trans - (cl/v) * central
    d/dt(trans)   =  ktr * depot - ktr * trans
    ## Concentration is calculated
    cp = central/v
    # And is assumed to follow proportional and additive error
    cp ~ prop(prop.err) + add(add.err)
  })
}

run006F <-
  nlmixr2(KAtr1_sexV, PKdata, est = "focei", foceiControl(print = 20))
run006F
save(run006F, file = "run006F.Rdata")
#load(file = "fitKAtr1_sexV_F.Rdata")


## binary categorical effects estimated on log scale can be back-transformed
THETA_SexV <- run006F$parFixedDf[4, ]
THETA_SexV
#          Parameter  Estimate        SE     %RSE Back-transformed    CI Lower  CI Upper BSV(CV%)
# Sex_V log Sex on v 0.4044965 0.2397131 59.26209        0.4044965 -0.06533252 0.8743256       NA
# Shrink(SD)%
# Sex_V          NA

#While it says 'backtransformed', this is not actually the case


## this results in fold-change estimates
exp(c(
  THETA_SexV$Estimate,
  THETA_SexV$"CI Lower",
  THETA_SexV$"CI Upper"))
#    Sex_V                     
#1.4985478 0.9367559 2.3972580 

## that can be translated in percentage change with proper confidence intervals
paste0(round(100 * (exp(c(
  THETA_SexV$Estimate,
  THETA_SexV$"CI Lower",
  THETA_SexV$"CI Upper")) - 1), 1), "%")
#[1] "49.9%"  "-6.3%"  "139.7%"



#################################################################################
##                                                                             ##
## Implement allometric covariates                                             ##
##                                                                             ##
#################################################################################

## Code is most stable if transformations are carried out in the data file
## instead of in the model code, especially for SAEM
## Using standard R syntax:
## PKdata$logWT70 <- log(PKdata$WT/70)

## Or using data.table syntax
PKdata[,logWT70:=log(WT/70)]

## One compartment transit model with allometric scaling on WT

One.comp.transit.allo <- function() {
  ini({
    # Where initial conditions/variables are specified
    lktr <- log(1.15)  #log k transit (/h)
    lcl  <- log(0.15)  #log Cl (L/hr)
    lv   <- log(7)     #log V (L)
    ALLC <- fix(0.75)  #allometric exponent cl
    ALLV <- fix(1.00)  #allometric exponent v
    prop.err <- 0.15   #proportional error (SD/mean)
    add.err <- 0.6     #additive error (mg/L)
    eta.ktr ~ 0.5
    eta.cl ~ 0.1
    eta.v ~ 0.1
  })
  model({
    #Allometric scaling on weight
    cl <- exp(lcl + eta.cl + ALLC * logWT70)
    v  <- exp(lv + eta.v + ALLV * logWT70)
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

run007F <-
  nlmixr(One.comp.transit.allo,
         PKdata,
         est = "focei",
         foceiControl(print = 5))
run007F


save(run007F, file = "run007F.Rdata")
#load(file = "run007F.Rdata")

## do you get a significant drop in OFV by including allometric weight?
load(file="run003F.Rdata")
run007FF$OBJF-run003F$OBJF
#[1] -29.98403

#################################################################################
##                                                                             ##
##  Hands-on assignments: nlmixr2 model development                             ##
##  Run the allometric model without fixing the exponents                      ##
##  Do you get a better fit?                                                   ##
##                                                                             ##
##  fitPK001$OBJF-fitPK002$OBJF                                                ##
##                                                                             ##
#################################################################################

#...

#################################################################################
##                                                                             ##
#################################################################################


