## load the required libraries
library(nlmixr2)
library(babelmixr2)
library(nonmem2rx)

library(lattice)


# babelmixr2 allows running nlmixr2 models using NONMEM as engine, 
# can provide starting values if the data set contains rich enough profiles to allow NCA analysis (using package PKNCA), 
# can keep track of units and much more, but that is material for another course!


#######################################################################################
#
# The nonmem2rx library allows you to read in NONMEM output to 
# create input for rxode2 or nlmixr2
#
#######################################################################################

# Read in NONMEM results -in this case from run005- into an rxode2 object
# Models read in correctly if Uppsala-style residual errors are used 
# with SIGMA fixed to 1 and with only a single SIGMA
# Example with combined additive and proportional error:

mod5 <- nonmem2rx("run005.ctl",lst=".res", save=FALSE)

# NONMEM results are read in along with the model that has been converted 
# to rxode2/nlmixr2 syntax with updated starting values
print(mod5)

# to check if the model has been implemented correctly, NONMEM PRED, IPRED, and IWRES are 
# compared with the rxode2 generated values after conversion
plot(mod5)

# the object can be used to simulate new regimens
# in this case four doses 12 hours apart of 100 mg using 100 subjects
ev <- et(amt=100, ii=12, until=48) %>%
  et(seq(0, 48, by=0.1)) %>%
  et(id=1:100)

s1 <- rxSolve(mod5, ev)
xyplot(center~time,data= s1, groups=id,type='l',xlab="Time (h)",ylab="Warfarin (mg/L)",
       scales = list(x=list(at=seq(0,60,6))))


# Or simulate 100 studies with these NONMEM parameters to capture uncertainty in theta estimates
s <- rxSolve(mod5, ev, nStud=100)
# calculate confidence interval estimates using 100*100 individual curves
sci <- confint(s, parm=c("gut","center"))
plot(sci)


# to continue with nlmixr2 functionality, translate the rxode2 object to an nlmixr2 object
run005NM <- as.nlmixr2(mod5)

# now you can generate individual concentration-time profiles for both PRED and IPRED
indivs5 <- augPred(run005NM)

## default plot
plot(indivs5)

## nicer plot with lattice

## specify array of colours for curves
nlmixCOLS <- c("#28466A","#8DB6CD","#B40000")  

xyplot(
  values ~ time/24 | id,
  data = indivs5,
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


# Create VPCs for the NONMEM model that has been translated to nlmixr2

#Find unique nominal sampling times
bin_mids <- sort(unique(round(as.numeric(run005NM$origData$TIME),3)))
#Put the edges in the middle of the nominal time points
bin_edges <- bin_mids - c(0, diff(bin_mids) / 2) 

vpc5NM<-vpcPlot( 
  run005NM,
  bins=bin_edges,
  n = 500,                       #number of replicates simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE,
              obs_median=TRUE,
              sim_median=TRUE,
              sim_median_ci=TRUE,
              obs_ci=TRUE,
              pi=TRUE),
  xlab = "Time (h)",             #x-axis label
  ylab = "Warfarin (mg/L)", #y-axis label
  title= "VPC for NONMEM run005"
)
vpc5NM


#results for the NONMEM-translated model
run005NM
#          Parameter   Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#theta1       log ka -0.569  0.246 43.3 -0.569 (-1.05, -0.0866)                     
#theta2       log cl  -2.01 0.0538 2.68     -2.01 (-2.11, -1.9)                     
#theta3        log v   2.06  0.046 2.23       2.06 (1.97, 2.15)                     
#EPS_prop   EPS_prop  0.117                               0.117                     
#EPS_pkadd EPS_pkadd  0.722                               0.722                     
#ETA_ka                                                             78.0      45.8% 
#ETA_cl                                                             28.7      3.74% 
#ETA_v                                                              22.2      13.4% 



run005F <-
  nlmixr2(run005NM,                 #the model object
          data=run005NM$origData,   #the data inside the model object
          est = "focei",            #the estimation algorithm (FOCEi)
          #FOCEi options:
          foceiControl(print = 5))  #only print every 5th estimation step




#nlmixr2 results after running the NONMEM-translated model
run005F
#           Parameter   Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
# theta1       log ka -0.569  0.344 60.4   -0.569 (-1.24, 0.105)                     
# theta2       log cl  -2.01 0.0776 3.87    -2.01 (-2.16, -1.86)                     
# theta3        log v   2.06 0.0669 3.25       2.06 (1.93, 2.19)                     
# EPS_prop   EPS_prop  0.117                               0.117                     
# EPS_pkadd EPS_pkadd  0.722                               0.722                     
# ETA_ka                                                             78.0      45.8% 
# ETA_cl                                                             28.7      3.74% 
# ETA_v                                                              22.2      13.4%                                                    
