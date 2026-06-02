################################################################################
##  Hands-on Session : Integrating nlmixr2autoinit with nlmixr2
##
##  nlmixr2autoinit Tutorial  |  UCL Pharmacometrics Group
##
##  Objectives:
##    A. Run getPPKinits() on warfarin PK data - explore nlmixr2autoinit outputs
##    B. Fit a 1-compartment oral model (SAEM) using nlmir2autoinit estimates
##    C. Run getPPKinits() on your own data and fit a 1-compartment model (optional)
################################################################################

## ── Setup ─────────────────────────────────────────────────────────────────────
# install.packages("nlmixr2autoinit", dependencies = TRUE)
library(nlmixr2)
library(nlmixr2autoinit)

################################################################################
## PART A - getPPKinits() on phenobarbital data (IV)
################################################################################

## ── Prepare warfarin PK data ──────────────────────────────────────────────────
# Read warfarin from nlmixr2data and keep concentration observations only
warf_pk <- warfarin[warfarin$dvid == "cp", ]
warf_pk$dvid <- NULL
rownames(warf_pk) <- NULL
warf_pk$sex <- ifelse(warf_pk$sex == "male", 1, 0)
# nlmixr2autoinit requires upper-case column names
colnames(warf_pk) <- toupper(colnames(warf_pk))

## ── Quick data plot ───────────────────────────────────────────────────────────
library(ggplot2)
obs_only <- warf_pk[warf_pk$EVID == 0, ]

ggplot(warf_pk, aes(x = TIME, y = DV, group = ID)) +
  geom_line(alpha = 0.3,size=0.8) + 
  # geom_point(size = 1) +
  scale_x_continuous("Time (h)") +
  scale_y_continuous("Warfarin concentration (mg/L)") +
  theme_bw() + theme(legend.position = "none") +
  theme(
    axis.text        = element_text(color = "#1a1a1a", size = 12),
    axis.title       = element_text(color = "#1a1a1a", size = 15))

################################################################################
## PART A - Run getPPKinits() on warfarin data
################################################################################

# Note: for an oral dataset, provide a cmt column to indicate the depot
# compartment; otherwise nlmixr2autoinit treats dosing as IV.
warf_pk$CMT<-2
warf_pk[warf_pk$EVID==1,]$CMT<-1

inits_warf <- getPPKinits(dat = warf_pk)

## A.1  Data information detected
inits_warf$Datainfo
# Expected: Dose Route = oral, Dose Type = first_dose, Subjects = 32

## A.2  Recommended initial estimates
inits_warf$Recommended_initial_estimates
# Ka, CL, Vd; also 2-/3-cmt and Vmax/Km results from parameter sweeps

## A.3  Method comparison for Ka / CL / Vd
inits_warf$Run.history$base.out
# metrics.rank = 1 = best; compare rRMSE2 across method combinations
# ?metrics.  # read metrics details

# Michaelis-Menten (Vmax/Km) parameter sweeping results
inits_warf$Run.history$sim.vmax.km

## A.4  Multi-compartment parameter sweep results
# 2-compartment parameter sweeping results (Q, V2/Vp estimates)
inits_warf$Run.history$sim.2cmpt   # inspect column names for Q / Vp / V2
# 3-compartment parameter sweeping results
inits_warf$Run.history$sim.3cmpt
## A.5  Residual error estimation
inits_warf$Run.history$sigma.out

################################################################################
## PART B - Fit 1-compartment oral model (SAEM)
################################################################################

## B.1  Extract initial estimates from nlmixr2autoinit output
est <- inits_warf$Recommended_initial_estimates

Ka       <- as.numeric(est$Values[est$Parameters == "Ka"])
CL       <- as.numeric(est$Values[est$Parameters == "CL"])
Vd       <- as.numeric(est$Values[est$Parameters == "Vd"])
add_err  <- as.numeric(est$Values[est$Parameters == "Sigma additive"])
prop_err <- as.numeric(est$Values[est$Parameters == "Sigma proportional"])

## B.2  Define 1-compartment oral model using nlmixr2autoinit estimates
warf_1cmt<- function() {
  ini({
    lka  <- log(Ka)        # log Ka (1/h) 
    lcl  <- log(CL)        # log CL (L/h) 
    lv   <- log(Vd)        # log Vd (L)   
    # IIV - default values
    eta.ka ~ 0.1
    eta.cl ~ 0.1
    eta.v  ~ 0.1
    # Residual error 
    prop.err <- prop_err
    add.err  <- add_err
  })
  model({
    ka <- exp(lka + eta.ka)
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv  + eta.v)
    
    d/dt(depot)   = -ka * depot
    d/dt(central) =  ka * depot - (cl/v) * central
    
    cp = central / v
    cp ~ prop(prop.err) + add(add.err)
  })
}

## B.3  Check model structure before fitting
nlmixr2(warf_1cmt)

## B.4  Fit with SAEM
run001_warf <- nlmixr2(
  warf_1cmt,
  warf_pk,
  est     = "saem",
  control = saemControl(nBurn = 200, nEm = 300, print = 50, logLik = TRUE),
  table   = tableControl(cwres = TRUE)
)

print(run001_warf)
saveRDS(run001_warf, file = "run001_warf.RDS")

# Check results
knitr::kable(run001_warf$parFixed, digits = 3,
      caption = "1-cmt model estimates (SAEM)") %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)


################################################################################
## PART C (Optional) - Use your own data
##
##  Replace pheno_sd below with your own dataset.
##  Two model templates provided: IV and oral 1-compartment.
##  Adjust parameter names, IIV structure, and error model as needed.
################################################################################

## C.1  Load your data (pheno_sd used here as placeholder)
# mydata <- read.csv("your_data.csv")   # <- replace with your file
mydata <- pheno_sd

## C.2  Run getPPKinits()
inits_my <- getPPKinits(dat = mydata)
inits_my$Datainfo
inits_my$Recommended_initial_estimates

## C.3  Extract estimates
est_my  <- inits_my$Recommended_initial_estimates
CL_my   <- as.numeric(est_my$Values[est_my$Parameters == "CL"])
Vd_my   <- as.numeric(est_my$Values[est_my$Parameters == "Vd"])
add_my  <- as.numeric(est_my$Values[est_my$Parameters == "Sigma additive"])
prop_my <- as.numeric(est_my$Values[est_my$Parameters == "Sigma proportional"])

## ── Template 1: IV 1-compartment ──────────────────────────────────────────
# Use if your data has IV dosing (bolus or infusion). No Ka needed.

my_1cmt_iv <- function() {
  ini({
    lcl <- log(CL_my)      # log CL (L/h) 
    lv  <- log(Vd_my)      # log Vd (L)  
    eta.cl ~ 0.1           # IIV on CL 
    eta.v  ~ 0.1           # IIV on Vd 
    prop.err <- prop_my
    add.err  <- add_my
  })
  model({
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv  + eta.v)
    d/dt(central) = -(cl/v) * central
    cp = central / v
    cp ~ prop(prop.err) + add(add.err)
  })
}

run_iv <- nlmixr2(
  my_1cmt_iv, mydata, est = "saem",
  control = saemControl(nBurn = 200, nEm = 300, print = 50, logLik = TRUE),
  table   = tableControl(cwres = TRUE)
)

saveRDS(run_iv, file = "run_my_iv.RDS")

knitr::kable(run_iv$parFixed, digits = 3,
             caption = "1-cmt model estimates (SAEM)") %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)

## ── Template 2: Oral 1-compartment ────────────────────────────────────────
# Use if your data has oral/extravascular dosing.
# Requires Ka; ensure your dataset has a CMT column for the depot compartment.

Ka_my <- as.numeric(est_my$Values[est_my$Parameters == "Ka"])  # oral only

my_1cmt_oral <- function() {
  ini({
    lka  <- log(Ka_my)     # log Ka (1/h) 
    lcl  <- log(CL_my)     # log CL (L/h) 
    lv   <- log(Vd_my)     # log Vd (L)  
    eta.ka ~ 0.1           # IIV on Ka 
    eta.cl ~ 0.1           # IIV on CL
    eta.v  ~ 0.1           # IIV on Vd
    prop.err <- prop_my
    add.err  <- add_my
  })
  model({
    ka <- exp(lka + eta.ka)
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv  + eta.v)
    d/dt(depot)   = -ka * depot
    d/dt(central) =  ka * depot - (cl/v) * central
    cp = central / v
    cp ~ prop(prop.err) + add(add.err)
  })
}

run_oral <- nlmixr2(
  my_1cmt_oral, mydata, est = "saem",
  control = saemControl(nBurn = 200, nEm = 300, print = 50, logLik = TRUE),
  table   = tableControl(cwres = TRUE)
)
print(run_oral)
saveRDS(run_oral, file = "run_my_oral.RDS")

knitr::kable(run_oral$parFixed, digits = 3,
             caption = "1-cmt model estimates (SAEM)") %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)

