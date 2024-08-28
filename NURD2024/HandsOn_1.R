## load the required libraries
library(nlmixr2)
library(lattice)
## (or use your favourite graphics package)


#################################################################################
##                                                                             ##
## rxode2 simulation Part 1                                                    ##
##                                                                             ##
## Simulation of a single (warfarin concentration) curve with a single dose    ##
##                                                                             ##
#################################################################################

## set up the system of differential equations (ODEs)
odeKA1 <- "
 d/dt(depot)   = -ka*depot;
 d/dt(central) =  ka*depot-(cl/v)*central;
 C1=central/v;
"

## compile the model
modKA1 <- rxode2(model = odeKA1)

## provide the parameter values to be simulated:
Params <- c(
 ka = log(2) / 0.5, # 1/h (aborption half-life of 30 minutes)
 cl = 0.135,        # L/h
 v = 8              # L
)

## create an empty event table that stores both dosing and sampling information :
ev <- eventTable()

## add a dose to the event table:
ev$add.dosing(dose = 500) #mg
              
## add time points to the event table where concentrations will be simulated; these actions are cumulative
ev$add.sampling(seq(0, 120, 0.1))

## Then solve the system
##
## The output from rxSolve is a solved RxODE object,
##  but by making it a data.frame only the simulated values are kept:
Res<-data.frame(rxSolve(modKA1,Params,ev))
#Alternative code:
Res<-data.frame(modKA1$run(Params,ev))


## then plot the simulated outcomes in the compartments:
## the amounts in the depot compartment
xyplot(depot~time,data=Res,type='l',lwd=2)
## the concentrations in the central compartment
xyplot(C1~time,data=Res,type='l',lwd=2)


## Extend the eventTable by adding three infusions to the central compartment
## Remember: updates to the eventTable are cumulative

ev$add.dosing(
  dose = 250,           #mg
  nbr.doses = 3,        #add three doses
  dosing.to = 2,        #add them to the second ODE in the model (=central)
  dosing.interval = 12, #h; set the doses 12 hours apart
  rate = 125,           #mg/h; infuse at a rate of 125 mg/h, resulting in 2-hour infusions
  start.time = 36       #h; have the three doses start at 36h
)

Res<-data.frame(rxSolve(modKA1,Params,ev))

## only the first dose goes into the depot compartment, so nothing changes here:
xyplot(depot~time,data=Res,type='l',lwd=2)

## the concentrations in the central compartment now reflect the three additional infusions: in the compartments:
p1<-xyplot(C1~time,data=Res,type='l',lwd=2)
p1



#################################################################################
##                                                                             ##
##  Hands-on assignments: rxode2 simulation                                    ##
##  -Simulate the effect of changing parameters                                ##
##  -Simulate the effect of changing doses                                     ##
##                                                                             ##
#################################################################################

#...

#################################################################################
##                                                                             ##
#################################################################################





#################################################################################
##                                                                             ##
## Extending the model with a transit compartment                              ##
##                                                                             ##
#################################################################################

## This requires an update to the ODE equation for the central compartment
##  because amounts now come from trans instead of depot:


odeKA1trans <- "
 d/dt(depot)   = -ka*depot;
 d/dt(central) =  ktr*trans-(cl/v)*central;
 d/dt(trans)   =  ka*depot-ktr*trans; 
 C1=central/v;
"

## compile the model
modKA1trans <- RxODE(model = odeKA1trans)

## provide the extra ktr parameter:
Params2 <- c(
  ka = log(2) / 0.5, # 1/h (aborption half-life of 30 minutes)
  cl = 0.135,        # L/h
  v = 8,             # L
  ktr = log(2) / 5   # 1/h (transit half-life of 5 hours)
)

## the eventTable does not have to change

Res2<-data.frame(rxSolve(modKA1trans,Params2,ev))

## Examine the results:
## Only the first dose goes into the depot compartment, so nothing changes here:
xyplot(depot~time,data=Res2,type='l',lwd=2)

## Only the first dose gets transferred into the transit compartment:
xyplot(trans~time,data=Res2,type='l',lwd=2)

## The concentration in the central compartment then becomes:
p2<-xyplot(C1~time,data=Res2,type='l',lwd=2)
p2


print(p1,position=c(0,0.5,1,1),more=TRUE)
print(p2,position=c(0,0,1,0.5),more=FALSE)






#################################################################################
##                                                                             ##
##  Alternatively you can provide a NONMEM-style data set                      ##
##  instead of an eventTable                                                   ##
##                                                                             ##
#################################################################################


library(data.table)
evNM<-data.table(time=seq(0, 120, 0.1))
evNM[,amt:=ifelse(time==0,500,NA)]
evNM[,cmt:=ifelse(time==0,1,2)]


evNM[,amt:=ifelse(time==36,250,amt)]
evNM[,rate:=ifelse(time==36,amt/2,NA)]
evNM[,ii:=ifelse(time==36,12,NA)]
evNM[,addl:=ifelse(time==36,2,NA)]

Res2NM<-data.frame(rxSolve(modKA1trans,Params2,evNM))

## Examine the results:
## Only the first dose goes into the depot compartment, so nothing changes here:
#xyplot(depot~time,data=Res2NM,type='l',lwd=2)

## Only the first dose gets transferred into the transit compartment:
#xyplot(trans~time,data=Res2NM,type='l',lwd=2)

## The concentration in the central compartment then becomes:
p2NM<-xyplot(C1~time,data=Res2NM,type='l',lwd=2)
p2NM

#Compare the two simulated curves
print(p2,position=c(0,0.5,1,1),more=TRUE)
print(p2NM,position=c(0,0,1,0.5),more=FALSE)







#################################################################################
##                                                                             ##
##  Hands-on bonus assignment: ODE update                                      ##
##  -Change the system of ODEs and examine the results                         ##
##  For example: extend the original model with five transit compartments      ##
##  and use 4 bolus doses in the 1st compartment                               ##
##  -Add an effect compartment and simulate an effect                          ##
##                                                                             ##
#################################################################################

#...

#################################################################################
##                                                                             ##
#################################################################################

