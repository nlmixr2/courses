####################################################################################
   ## Project   : Generic code to simulate exposures by weight and age
   ##              for different dosing algorithms 
   ## Purpose   : Simulations and graphs
   ## Client    : Occams Cooperatie UA
   ## Reference : Nhanes101
##### ------------------------------------------------------------------------------
   ## Author: 
   ## Rik Schoemaker, PhD
   ## Occams
   ## rik.schoemaker@occams.com
   ## +31 20 441 6410
##### ------------------------------------------------------------------------------
   ## Copyright (C) 2024 Occams
   ## The code in this script is free for anyone to use.
#####################################################################################

## remove all pre-existing information
rm(list=ls())

## load required libraries
library(rxode2)
library(ggplot2)
library(patchwork)
library(mvnfast)
library(ggforce)
library(nonmem2rx)
library(data.table)
library(lattice)




#######################################################################################
##
## Paediatric PopPK model run625
##
#######################################################################################

## Read in EBEs for kids in the analyses

EBEsP<-fread("run625.csv")
## create a variable of age in months
EBEsP[,Age:=AGE*12]

## Four dosing scenarios to examine: 2mg/kg bid, 2.5mg/kg bid and 3 mg/kg bid all with a maximimum of 100mg bid
## followed by a scenario split in weight groups:
## <10kg: 3 mg/kg bid, 10-<20kg: 2.5 mg/kg bid, >=20kg: 2 mg/kg bid with a maximum of 100 mg bid

EBEsP[,AMT2mgkg:=2*WT]
EBEsP[,AMT2mgkg:=ifelse(AMT2mgkg>100,100,AMT2mgkg)]
EBEsP[,AMT25mgkg:=2.5*WT]
EBEsP[,AMT25mgkg:=ifelse(AMT25mgkg>100,100,AMT25mgkg)]
EBEsP[,AMT3mgkg:=3*WT]
EBEsP[,AMT3mgkg:=ifelse(AMT3mgkg>100,100,AMT3mgkg)]
EBEsP[,AMT3252:=ifelse(WT<10,AMT3mgkg,
                          ifelse(WT<20,AMT25mgkg,AMT2mgkg))]

## Create a sampling schedule for bid dosing; first 2 doses and and then a total of 24 doses to reach Steady state
## Time interval for sampling:
StepSize<-0.1
## Times post-dose (over 12 hours because bid dosing:
timeS<-c(0.001,seq(StepSize,12,StepSize))
## Full sampling profile with profiles starting at 0h, 12h and 23*12h 
timeZ<-c(0,timeS,12+timeS,276+timeS)

## Use data.table syntax to combine the ID values in EBEsP with the time schedule to create a NONMEM type data structure
## with all ID and time combinations:
SIMP<-EBEsP[,.(TIME=timeZ),.(ID)]
## then merge the structure with the EBEs and soing information:
SIMP<-merge(SIMP,EBEsP,by="ID")
## set up AMT records at time 0 
SIMP[,AMT2:=ifelse(TIME ==0,AMT2mgkg,0)]
SIMP[,AMT25:=ifelse(TIME ==0,AMT25mgkg,0)]
SIMP[,AMT3:=ifelse(TIME ==0,AMT3mgkg,0)]
SIMP[,AMT3252:=ifelse(TIME ==0,AMT3252,0)]
## make all records that are not dosing (at time zero) into 'other type' event records:
SIMP[,EVID:=ifelse(TIME== 0,1,2)]
## set the bid dosing interval to 12:
SIMP[,II:=ifelse(TIME== 0,12,0)]
## create 24 additional doses using the ADDL data item:
SIMP[,ADDL:=ifelse(TIME== 0,24,0)]
## this NONMEM type data structure withe EBEs can now be used to simulate individual profiles

## Define the ODE for the NONMEM model
ode625<-"
 k10 <- CL/VC
 d/dt(gut) = -KA*gut
 d/dt(centr) = KA*gut-k10*centr
 ipred <- centr / VC
"
#Compile the model
mode625 <- rxode2(model=ode625)

## Simulate 2mg bid with 100mg bid maximum by setting AMT to the correct value 
SIMP[,AMT:=AMT2]
## Then simulate all profiles by combining the model and the data structure using rxSolve
Sim2P <- as.data.table(rxSolve(mode625, SIMP))

## Simulate 2.5mg bid with 100mg bid maximum
SIMP[,AMT:=AMT25]
Sim25P <- as.data.table(rxSolve(mode625, SIMP))

## Simulate 3mg bid with 100mg bid maximum
SIMP[,AMT:=AMT3]
Sim3P <- as.data.table(rxSolve(mode625, SIMP))

## Simulate proposed posology
SIMP[,AMT:=AMT3252]
Sim3252P <- as.data.table(rxSolve(mode625, SIMP))



## Calculate summary measures on the data sets and merge them with the EBEs and dosing information:
## Adjust this to generate the required exposure measures for each specific case
CalcSum<-function(Data=Sim2P,Demo=EBEsP){
  ## Cav at steady state:
  SumCavss<-Data[time>276,.(Cavss=mean(ipred)),.(id)]
  ## Cmax first dose:
  SumCmax1<-Data[time<12,.(Cmax1=max(ipred)),.(id)]
  ## trough at 24h:
  SumTrough24<-Data[time==24,.(Ctrough24=ipred),.(id)]
  SumP<-merge(SumCavss,SumCmax1,by='id')
  SumP<-merge(SumP,SumTrough24,by='id')
  SumP<-merge(SumP,Demo,by='id')
  SumP
}


## Summary measures for paediatric EBEs:

## rxode2 has the peculiarity to change all variable names to lower case; 
## create a lower case version of the ID variable for merging
EBEsP[,id:=ID]

SumEBEsP2<-CalcSum(Sim2P,EBEsP)
SumEBEsP25<-CalcSum(Sim25P,EBEsP)
SumEBEsP3<-CalcSum(Sim3P,EBEsP)
SumEBEsP3252<-CalcSum(Sim3252P,EBEsP)








#######################################################################################
##
## Simulate paediatric subjects using the Nhanes database and PopPK model run625
##
#######################################################################################

## Simulations of exposures where demographic variables are sampled from the Nhanes database
## Usually only weight and age are used, but a model with lean body weight can be driven by this dataset as well

Nhanes<- fread("NhanesDemo2006.csv")

## Paediatric subjects:
NhanesP<-Nhanes[AGE<18]

## Reduce the dataset by sampling from 4 month age bins for 0-4 years, and 1 year age bins for 4-18 years with 200 subjects per bin
## This example focuses on 0-48 months only and the extra bins from 4-18 years are for illustrative purposes only
NhanesP[,AgeCat:=cut(AGE,breaks=c(seq(0,4,0.25),seq(5,18),seq(28,88,10)),right=FALSE,include.lowest=TRUE)]

## Make a variable of age in months
NhanesP[,Age:=AGE*12]

set.seed(183476)
NhanesP<-NhanesP[,.SD[sample(.N,200,replace=TRUE)],.(AgeCat)]

## Create ID variable for simulation
NhanesP[,ID:=1:.N]

## Add AED covariate by sampling from a uniform distribution with the observed AED fraction in kids
## Determine the fraction of subjects taking AED 
FAEDP<-unlist(EBEsP[,.(mean(as.numeric(AED==1)))])
FAEDP
#1: 0.3957447

set.seed(183476)
NhanesP[,RUNIF:=runif(.N)]
NhanesP[,AED:=as.numeric(RUNIF<=FAEDP)]
NhanesP[,.(mean(AED))]
#1:  0.39
## so AED fraction is implemented correctly

#Potential dosing schedules:
NhanesP[,AMT2mgkg:=2*WT]
NhanesP[,AMT2mgkg:=ifelse(AMT2mgkg>100,100,AMT2mgkg)]
NhanesP[,AMT25mgkg:=2.5*WT]
NhanesP[,AMT25mgkg:=ifelse(AMT25mgkg>100,100,AMT25mgkg)]
NhanesP[,AMT3mgkg:=3*WT]
NhanesP[,AMT3mgkg:=ifelse(AMT3mgkg>100,100,AMT3mgkg)]
NhanesP[,AMT3252:=ifelse(WT<10,AMT3mgkg,
                         ifelse(WT<20,AMT25mgkg,AMT2mgkg))]




#######################################################################################
##
## Two ways are described to generate individual simulated profiles
## The first simulates individual parameters from the NONMEM estimates 
## and puts them in the simulation file
##
#######################################################################################



## extract the parameter estimates from the run625.res output

# TVCL   = THETA(2)
# TVVC   = THETA(3)
# TVKA   = THETA(4)
# ALLOCL = THETA(5)
# ALLOVC = THETA(6)
# AEDCL  = THETA(7)
# 
# CL    = EXP(LOG(TVCL) + ALLOCL*(LOG(WT/70)) + AEDCL*AED + ETA(1)) ;Clearance (L/hr)
# VC    = EXP(LOG(TVVC) + ALLOVC*(LOG(WT/70)) + ETA(2)) ;Volume (L)
# KA    = EXP(LOG(TVKA) + ETA(3)) ;Ka (1/hr) 


THETAP<-c(3.11E-01,  4.19E+00,  6.49E+01,  2.10E+00,  7.50E-01,  1.00E+00,  2.36E-01)
  
OMEGAP<-matrix(
  c( 5.31E-02,  0.00E+00, 0.00E+00,
     0.00E+00,  1.88E-02, 0.00E+00,
     0.00E+00,  0.00E+00, 2.93E-09),ncol=3)

## simulate individual ETAs from the NONMEM-estimated Omega matrix and name the columns ET1, ETA2 etc
set.seed(2093875)
ETAP <-
  as.data.table(rmvn(n = dim(NhanesP)[1], rep(0, dim(OMEGAP)[1]), OMEGAP)) # Sample from covariance matrix
names(ETAP)<-paste0("ETA",seq(1,dim(ETAP)[2]))

## bind the ETAs to the Nhanes subjects
NhanesPX<-cbind(NhanesP,ETAP)

## simulate PK parameters using the formulas in run625
NhanesPX[,CL:=exp(log(THETAP[2]) + THETAP[5]*(log(WT/70)) + THETAP[7]*AED + ETA1)]
NhanesPX[,VC:=exp(log(THETAP[3]) + THETAP[6]*(log(WT/70))+ETA2)]
NhanesPX[,KA:=exp(log(THETAP[4])+ETA3)]

## Make a NONMEM-like data structure with observation times and dosing records and covariates
## similar to what was done above for the kids in the trials
NhanesPSX<-NhanesPX[,.(TIME=timeZ),.(ID)]
NhanesPSX<-merge(NhanesPSX,NhanesPX,by="ID")

## make the dose records
NhanesPSX[,AMT2:=ifelse(TIME ==0,AMT2mgkg,0)]
NhanesPSX[,AMT25:=ifelse(TIME ==0,AMT25mgkg,0)]
NhanesPSX[,AMT3:=ifelse(TIME ==0,AMT3mgkg,0)]
NhanesPSX[,AMT3252:=ifelse(TIME ==0,AMT3252,0)]
NhanesPSX[,EVID:=ifelse(TIME== 0,1,2)]
NhanesPSX[,II:=ifelse(TIME== 0,12,0)]
NhanesPSX[,ADDL:=ifelse(TIME== 0,24,0)]


## Simulate 2mg/kg bid with 100mg bid maximum using the created data file with EBEs and the ode defined above (mode625)
NhanesPSX[,AMT:=AMT2]
PredP2X<-as.data.table(rxSolve(mode625,NhanesPSX))

## Simulate 2.5mg/kg bid with 100mg bid maximum using the created data file with EBEs and the ode defined above (mode625)
NhanesPSX[,AMT:=AMT25]
PredP25X<-as.data.table(rxSolve(mode625,NhanesPSX))

## Simulate 3mg/kg bid with 100mg bid maximum using the created data file with EBEs and the ode defined above (mode625)
NhanesPSX[,AMT:=AMT3]
PredP3X<-as.data.table(rxSolve(mode625,NhanesPSX))

## Simulate proposed posology using the created data file with EBEs and the ode defined above (mode625)
NhanesPSX[,AMT:=AMT3252]
PredP3252X<-as.data.table(rxSolve(mode625,NhanesPSX))

## create a lower case id/ID variable
NhanesP[,id:=ID]

## Calculate summary measures
SumP2X<-CalcSum(PredP2X,NhanesP)
SumP25X<-CalcSum(PredP25X,NhanesP)
SumP3X<-CalcSum(PredP3X,NhanesP)
SumP3252X<-CalcSum(PredP3252X,NhanesP)






#######################################################################################
##
## Simulate adult subjects using the Nhanes database and PopPK model run025
##
#######################################################################################

## Simulations of exposures where demographic variables are sampled from the Nhanes database
## Usually only weight and age are used, but a model with lean body weight can be driven by this dataset as well

Nhanes<- fread("NhanesDemo2006.csv")

## Adult subjects:
NhanesA<-Nhanes[AGE>=18]

## Reduce the dataset by sampling 1000 subjects

set.seed(183476)
NhanesA<-NhanesA[,.SD[sample(.N,1000,replace=TRUE)],]

## Create ID variable for simulation
NhanesA[,ID:=1:.N]

## Add AED covariate by sampling from a uniform distribution with the observed AED fraction in Adults
## Determine the fraction of subjects taking AED 

EBEsA<-fread("run025.csv")
FAEDA<-unlist(EBEsA[,.(mean(as.numeric(AED==1)))])
FAEDA
# 0.5272331

set.seed(183476)
NhanesA[,RUNIF:=runif(.N)]
NhanesA[,AED:=as.numeric(RUNIF<=FAEDA)]
NhanesA[,.(mean(AED))]
#1: 0.515
## so AED fraction is implemented correctly

#flat dosing 100mg bid:
NhanesA[,AMT100:=100]

#######################################################################################
##
## Two ways are described to generate individual simulated profiles
## The first simulates individual parameters from the NONMEM estimates 
## and puts them in the simulation file
##
#######################################################################################

## extract the parameter estimates from the run025.res output

# TVCL   = THETA(2)
# TVVC   = THETA(3)
# TVKA   = THETA(4)
# ALLOCL = THETA(5)
# ALLOVC = THETA(6)
# AEDCL  = THETA(7)
# 
# CL    = EXP(LOG(TVCL) + ALLOCL*(LOG(WT/70)) + AEDCL*AED + ETA(1)) ;Clearance (L/hr)
# VC    = EXP(LOG(TVVC) + ALLOVC*(LOG(WT/70)) + ETA(2)) ;Volume (L)
# KA    = EXP(LOG(TVKA) + ETA(3)) ;Ka (1/hr) 


THETAA<-c(2.08E-01,  3.61E+00,  5.11E+01,  1.39E+00,  7.50E-01,  1.00E+00,  3.25E-01)

OMEGAA<-matrix(
  c( 6.90E-02,  0.00E+00, 0.00E+00,
     0.00E+00,  1.06E-01, 0.00E+00,
     0.00E+00,  0.00E+00, 7.93E-01),ncol=3)

## simulate individual ETAs from the NONMEM-estimated Omega matrix and name the columns ETA1, ETA2 etc
set.seed(2093875)
ETAA <-
  as.data.table(rmvn(n = dim(NhanesA)[1], rep(0, dim(OMEGAA)[1]), OMEGAA)) # Sample from covariance matrix
names(ETAA)<-paste0("ETA",seq(1,dim(ETAA)[2]))

## bind the ETAs to the Nhanes subjects
NhanesAX<-cbind(NhanesA,ETAA)

## simulate PK parameters using the formulas in run625
NhanesAX[,CL:=exp(log(THETAA[2]) + THETAA[5]*(log(WT/70)) + THETAA[7]*AED + ETA1)]
NhanesAX[,VC:=exp(log(THETAA[3]) + THETAA[6]*(log(WT/70))+ETA2)]
NhanesAX[,KA:=exp(log(THETAA[4])+ETA3)]

## Make a NONMEM-like data structure with observation times and dosing records and covariates
## similar to what was done above for the kids in the trials
NhanesASX<-NhanesAX[,.(TIME=timeZ),.(ID)]
NhanesASX<-merge(NhanesASX,NhanesAX,by="ID")

## make the dose records
NhanesASX[,AMT:=ifelse(TIME ==0,AMT100,0)]
NhanesASX[,EVID:=ifelse(TIME== 0,1,2)]
NhanesASX[,II:=ifelse(TIME== 0,12,0)]
NhanesASX[,ADDL:=ifelse(TIME== 0,24,0)]


## Simulate 100MG bid using the created data file with EBEs and the ode defined above (mode625);
## the ODE is the same for adults and children
PredAX<-as.data.table(rxSolve(mode625,NhanesASX))

## create a lower case id/ID variable
NhanesA[,id:=ID]

## Calculate summary measures
SumAX<-CalcSum(PredAX,NhanesA)



## simulated 90% of adult Cav values at 100 mg bid - 200 mg/day:
q5_640X<-quantile(SumAX$Cavss,0.05)
q95_640X<-quantile(SumAX$Cavss,0.95)


#######################################################################################
##
## Calculate the area containing 90% of the data using the function CalcArea()
##
#######################################################################################

CalcArea<-function(par="Cavss",by="Age",breaks=seq(0,48,4),data=SumP3252P){
  ## make a local copy
  Data<-copy(data)
  ## make sure that the data range does not extend beyond the breaks
  Data<-Data[get(by)>=breaks[1]&get(by)<=breaks[length(breaks)]]
  ## apply the breaks and create variables for the lower edge (BRKL) and the upper edge (BRKh) of the break interval
  Data[,BRKS:=cut(get(by),breaks=breaks,include.lowest=TRUE,right=TRUE)]
  Data[,BRKl:=cut(get(by),breaks=breaks,include.lowest=TRUE,right=TRUE,
                  labels=breaks[-length(breaks)])]
  Data[,BRKh:=cut(get(by),breaks=breaks,include.lowest=TRUE,right=TRUE,
                  labels=breaks[-1])]
  ## calculate summaries for the parameter (par) by the 'by' variable
  Sim<-Data[,.(Mean=mean(get(par)),
               Median=median(get(par)),
               Q2.5=quantile(get(par),probs=0.025),
               Q5=quantile(get(par),probs=0.05),
               Q25=quantile(get(par),probs=0.25),
               Q75=quantile(get(par),probs=0.75),
               Q95=quantile(get(par),probs=0.95),
               Q97.5=quantile(get(par),probs=0.975)),
            keyby=.(BRKS,BRKl,BRKh)]
  
  ## this creates a file where the low quantiles (5%) are set to run from the lower to the upper edge of the break, 
  ## and when you reach the end of the lower part of the area, it continues up to the high quantiles (95%) 
  ## and continues in the reverse direction, describing the area containing 90% of the data
  Sim$Area50<-Sim$Q25
  Sim$Area90<-Sim$Q5
  Sim$Area95<-Sim$Q2.5
  CssRev<-Sim[order(Sim$BRKS,decreasing=TRUE),]
  CssRev$Area50<-CssRev$Q75
  CssRev$Area90<-CssRev$Q95
  CssRev$Area95<-CssRev$Q97.5
  Css1<-Sim
  Css1$BY<-Css1$BRKl
  Css1$low<-1
  Css2<-Sim
  Css2$BY<-Css2$BRKh
  Css2$low<-0
  Css2<-rbind(Css1,Css2)
  Css2<-Css2[order(Css2$BY,Css2$low),]
  Css3<-CssRev
  Css3$BY<-Css3$BRKl
  Css3$low<-1
  Css4<-CssRev
  Css4$BY<-Css4$BRKh
  Css4$low<-0
  Css4<-rbind(Css3,Css4)
  Css4<-Css4[order(Css4$BY,Css4$low,decreasing=TRUE),]
  Css5<-rbind(Css2,Css4)
  Css5[,BY:=as.numeric(as.character(BY))]
  setnames(Css5,"BY",by)
  Css5
}



## specify the colours of the graph
AreaColour1<-"#d2edff"
AreaColour2<-"grey90"
AreaLineColour<-"black"
LineColour<-"dodgerblue"
MarkerColour<-"brown"

## calculate the area, in the case for the average steady state concentration over 0-48 months by weight
Cav3252mgWT<-CalcArea(par="Cavss",
                        by="WT",
                        breaks=c(3.4,6,8,10,12,14,16,18,20,28.1),
                        data=SumP3252X[Age<=48])

ps4WT<-xyplot(Cavss~WT,data=SumEBEsP3252[Age<=48],xlim=c(2.9,21.4),ylim=c(0.5,4.3),
              ylab=expression('Drug X C'[av]*' (mg/L)'),
              xlab="Weight (kg)",
              scales=list(alternating=1,y=list(at=seq(1,5,1)),x=list(at=seq(2,20,2))),
              main="Proposed posology",
              panel=function(x,y){
                tmp<-Cav3252mgWT
                lrect(-10,q5_640X,200,q95_640X,col=AreaColour2,border=AreaLineColour)
                panel.polygon(tmp$WT,tmp$Area90,col=AreaColour1,border=AreaLineColour)
                panel.xyplot(x,y,pch=16,cex=1.0,col=MarkerColour)
                panel.xyplot(tmp$WT,tmp$Median,type='l',col=LineColour,lwd=3)
                panel.abline(v=c(10,20))
              })
print(ps4WT)










#######################################################################################
##
## Two ways are described to generate individual simulated profiles
## The second simulates individual parameters by reading in the NONMEM output
## using the nonmem2rx package and then simulating from a data structure 
## that only contains covariates but not EBEs like in the first approach
##
#######################################################################################

## Read NONMEM paediatric model and parameter estimates to allow simulation using nonmem2rx
NM625 <- nonmem2rx("run625.ctl",lst=".res", save=FALSE)
## this allows an examination of the model 
## and setting of the population parameters as estimated by NONMEM:
print(NM625)

## Make a NONMEM-like data structure with observation times and dosing records and covariates
## similar to what was done above for the kids in the trials
NhanesPS<-NhanesP[,.(TIME=timeZ),.(ID)]
NhanesPS<-merge(NhanesPS,NhanesP,by="ID")

## make the dose records
NhanesPS[,AMT2:=ifelse(TIME ==0,AMT2mgkg,0)]
NhanesPS[,AMT25:=ifelse(TIME ==0,AMT25mgkg,0)]
NhanesPS[,AMT3:=ifelse(TIME ==0,AMT3mgkg,0)]
NhanesPS[,AMT3252:=ifelse(TIME ==0,AMT3252,0)]
NhanesPS[,EVID:=ifelse(TIME== 0,1,2)]
NhanesPS[,II:=ifelse(TIME== 0,12,0)]
NhanesPS[,ADDL:=ifelse(TIME== 0,24,0)]

## no individual PK parameters are required
## Simulate 2mg/kg bid with 100mg bid maximum using the created data file without EBEs and the model read using nonmem2rx (NM625)
NhanesPS[,AMT:=AMT2]
PredP2<-as.data.table(rxSolve(NM625,NhanesPS))
SumP2<-CalcSum(PredP2,NhanesP)

## Simulate 2.5mg/kg bid with 100mg bid maximum using the created data file without EBEs and the model read using nonmem2rx (NM625)
NhanesPS[,AMT:=AMT25]
PredP25<-as.data.table(rxSolve(NM625,NhanesPS))
SumP25<-CalcSum(PredP25,NhanesP)

## Simulate 3mg/kg bid with 100mg bid maximum using the created data file without EBEs and the model read using nonmem2rx (NM625)
NhanesPS[,AMT:=AMT3]
PredP3<-as.data.table(rxSolve(NM625,NhanesPS))
SumP3<-CalcSum(PredP3,NhanesP)

## Simulate proposed posology using the created data filewithout EBEs and the model read using nonmem2rx (NM625)
NhanesPS[,AMT:=AMT3252]
PredP3252<-as.data.table(rxSolve(NM625,NhanesPS))
SumP3252<-CalcSum(PredP3252,NhanesP)




#######################################################################################
##
## The same can be done with the adult data and model run025
##
#######################################################################################



## Read NONMEM adult model and parameter estimates to allow simulation using nonmem2rx
NM025 <- nonmem2rx("run025.ctl",lst=".res", save=FALSE)
## this allows an examination of the model 
## and setting of the population parameters as estimated by NONMEM:
print(NM025)

## Make a NONMEM-like data structure with observation times and dosing records and covariates
## similar to what was done above for the kids in the trials
NhanesAS<-NhanesA[,.(TIME=timeZ),.(ID)]
NhanesAS<-merge(NhanesAS,NhanesA,by="ID")

## make the dose records
NhanesAS[,AMT:=ifelse(TIME ==0,AMT100,0)]
NhanesAS[,EVID:=ifelse(TIME== 0,1,2)]
NhanesAS[,II:=ifelse(TIME== 0,12,0)]
NhanesAS[,ADDL:=ifelse(TIME== 0,24,0)]

## no individual PK parameters are required
## Simulate 100mg bid using the created data file without EBEs and the model read using nonmem2rx (NM025)
PredA<-as.data.table(rxSolve(NM025,NhanesAS))
SumA<-CalcSum(PredA,NhanesA)


#simulated 90% of adult Cav values at 100 mg bid - 200 mg/day:
q5_640<-quantile(SumA$Cavss,0.05)
q95_640<-quantile(SumA$Cavss,0.95)



Cav25mg<-CalcArea(par="Cavss",
                  by="Age",
                  breaks=seq(0,48,4),
                  data=SumP25[Age<=48])




AreaColour1<-"#d2edff"
AreaColour2<-"grey90"
AreaLineColour<-"black"
LineColour<-"dodgerblue"
MarkerColour<-"brown"
  
ps2<-xyplot(Cavss~Age,data=SumEBEsP25[Age<=48],xlim=c(-1,49),ylim=c(0.5,4.3),
            ylab=expression('Drug X C'[av]*' (mg/L)'),
            xlab="Age (months)",
            scales=list(alternating=1,x=list(at=seq(0,48,4)),y=list(at=seq(1,5,1))),
            main="2.5 mg/kg bid",
            panel=function(x,y){
              tmp<-Cav25mg
              lrect(-10,q5_640,200,q95_640,col=AreaColour2,border=AreaLineColour)
              panel.polygon(tmp$Age,tmp$Area90,col=AreaColour1,border=AreaLineColour)
              panel.xyplot(x,y,pch=16,cex=1,col=MarkerColour)
              panel.xyplot(tmp$Age,tmp$Median,type='l',col=LineColour,lwd=3)
            })
print(ps2)



#######################################################################################
##
## Alternatively, create the graph using ggplot2
##
#######################################################################################


AxisTextSize<-12
AxisTitleSize<-12
PlotTitleSize<-14

themeO<-theme(plot.margin=margin(0.2,0.2,0.2,0.2,"in"),
              panel.background=element_rect(fill = "white",colour="white"),
              axis.line = element_line(linewidth = 0.5, colour = "black"),
              axis.text=element_text(face='plain',size=AxisTextSize,colour='black'),
              axis.text.y=element_text(margin=margin(0,AxisTextSize/2,0,0,"pt")),
              axis.text.x=element_text(margin=margin(AxisTextSize/2,0,0,0,"pt")),
              axis.title.x=element_text(face='plain',margin=margin(AxisTitleSize,0,0,0,"pt"),size=AxisTitleSize),
              axis.title.y=element_text(face='plain',margin=margin(0,AxisTitleSize,0,0,"pt"),size=AxisTitleSize),
              axis.ticks.length=unit(AxisTextSize/2,'pt'),
              plot.title=element_text(face='plain',margin=margin(t = 0, r = 0, b = PlotTitleSize, l = 0, unit= "pt"),
                                      hjust=0.5,size=PlotTitleSize)
)


pg2<-ggplot()+
  geom_polygon(data=data.frame(x=c(-10,200,200,-10),y=c(q5_640,q5_640,q95_640,q95_640)),
               mapping=aes(x=x,y=y),
               fill=AreaColour2,colour=AreaLineColour)+
  geom_polygon(data=Cav25mg,mapping=aes(x=Age,y=Area90),fill=AreaColour1,colour=AreaLineColour)+
  geom_point(data=SumEBEsP25[Age<=48],mapping=aes(x=Age,y=Cavss),colour=MarkerColour,size=3,shape=16,stroke=0.8)+
  geom_path(data=Cav25mg,mapping=aes(x=Age,y=Median),colour=LineColour,linewidth=1.2)+
  scale_x_continuous(breaks=seq(0,48,4))+
  coord_cartesian(ylim=c(0.5,4.3),xlim=c(-1,49),expand=FALSE)+
  labs(y=expression('Drug X C'[av]*' (mg/L)'),
       x="Age (months)",
       title="2.5 mg/kg bid")+
  themeO

pg2




Cav3252mg<-CalcArea(par="Cavss",
                  by="Age",
                  breaks=seq(0,48,4),
                  data=SumP3252[Age<=48])


ps4<-xyplot(Cavss~Age,data=SumEBEsP3252[Age <=48],xlim=c(-1,49),ylim=c(0.5,4.3),
            ylab=expression('Drug X C'[av]*' (mg/L)'),
            xlab="Age (months)",
            scales=list(alternating=1,x=list(at=seq(0,48,4)),y=list(at=seq(1,5,1))),
            main="Proposed posology",
            panel=function(x,y){
              tmp<-Cav3252mg
              lrect(-10,q5_640,200,q95_640,col=AreaColour2,border=AreaLineColour)
              panel.polygon(tmp$Age,tmp$Area90,col=AreaColour1,border=AreaLineColour)
              panel.xyplot(x,y,pch=16,cex=1,col=MarkerColour)
              panel.xyplot(tmp$Age,tmp$Median,type='l',col=LineColour,lwd=3)
            })
print(ps4)


pg4<-ggplot()+
  geom_polygon(data=data.frame(x=c(-10,200,200,-10),y=c(q5_640,q5_640,q95_640,q95_640)),
               mapping=aes(x=x,y=y),
               fill=AreaColour2,colour=AreaLineColour)+
  geom_polygon(data=Cav3252mg,mapping=aes(x=Age,y=Area90),fill=AreaColour1,colour=AreaLineColour)+
  geom_point(data=SumEBEsP3252[Age<=48],mapping=aes(x=Age,y=Cavss),colour=MarkerColour,size=3,shape=16,stroke=0.8)+
  geom_path(data=Cav3252mg,mapping=aes(x=Age,y=Median),colour=LineColour,linewidth=1.2)+
  scale_x_continuous(breaks=seq(0,48,4))+
  coord_cartesian(ylim=c(0.5,4.3),xlim=c(-1,49),expand=FALSE)+
  labs(y=expression('Drug X C'[av]*' (mg/L)'),
       x="Age (months)",
       title="Proposed posology")+
  themeO

pg4

#Using the patchwork package
pg2+pg4




Cav25mgWT<-CalcArea(par="Cavss",
                    by="WT",
                    breaks=c(3.4,6,8,10,12,14,16,18,20,28.1),
                    data=SumP25[Age<=48])

ps2WT<-xyplot(Cavss~WT,data=SumEBEsP25[Age<=48],xlim=c(2.9,21.4),ylim=c(0.5,4.3),
              ylab=expression('Drug X C'[av]*' (mg/L)'),
              xlab="Weight (kg)",
              scales=list(alternating=1,y=list(at=seq(1,5,1)),x=list(at=seq(2,20,2))),
              main="2.5 mg/kg bid",
              panel=function(x,y){
                tmp<-Cav25mgWT
                lrect(-10,q5_640,200,q95_640,col=AreaColour2,border=AreaLineColour)
                panel.polygon(tmp$WT,tmp$Area90,col=AreaColour1,border=AreaLineColour)
                panel.xyplot(x,y,pch=16,cex=1.0,col=MarkerColour)
                panel.xyplot(tmp$WT,tmp$Median,type='l',col=LineColour,lwd=3)
              })
print(ps2WT)



pg2WT<-ggplot()+
  geom_polygon(data=data.frame(x=c(-10,200,200,-10),y=c(q5_640,q5_640,q95_640,q95_640)),mapping=aes(x=x,y=y),
               fill=AreaColour2,colour=AreaLineColour)+
  geom_polygon(data=Cav25mgWT,mapping=aes(x=WT,y=Area90),
               fill=AreaColour1,colour=AreaLineColour)+
  geom_point(data=SumEBEsP25[Age<=48],mapping=aes(x=WT,y=Cavss),
             colour=MarkerColour,size=3,shape=16,stroke=0.8)+
  geom_path(data=Cav25mgWT,mapping=aes(x=WT,y=Median),
            colour=LineColour,linewidth=1.2)+
  scale_x_continuous(breaks=seq(2,20,2))+
  coord_cartesian(ylim=c(0.5,4.3),xlim=c(2.9,21.4),expand=FALSE)+
  labs(y=expression('Drug X C'[av]*' (mg/L)'),
       x="Weight (kg)",
       title="2.5 mg/kg bid")+
  themeO

pg2WT


Cav3_25_2mgWT<-CalcArea(par="Cavss",
                        by="WT",
                        breaks=c(3.4,6,8,10,12,14,16,18,20,28.1),
                        data=SumP3252[Age<=48])

ps4WT<-xyplot(Cavss~WT,data=SumEBEsP3252[Age<=48],xlim=c(2.9,21.4),ylim=c(0.5,5.3),
              ylab=expression('Drug X C'[av]*' (mg/L)'),
              xlab="Weight (kg)",
              scales=list(alternating=1,y=list(at=seq(1,5,1)),x=list(at=seq(2,20,2))),
              main="Proposed posology",
              panel=function(x,y){
                tmp<-Cav3_25_2mgWT
                lrect(-10,q5_640,200,q95_640,col=AreaColour2,border=AreaLineColour)
                panel.polygon(tmp$WT,tmp$Area90,col=AreaColour1,border=AreaLineColour)
                panel.xyplot(x,y,pch=16,cex=1.0,col=MarkerColour)
                panel.xyplot(tmp$WT,tmp$Median,type='l',col=LineColour,lwd=3)
                panel.abline(v=c(10,20))
              })
print(ps4WT)

pg4WT<-ggplot()+
  geom_polygon(data=data.frame(x=c(-10,200,200,-10),y=c(q5_640,q5_640,q95_640,q95_640)),mapping=aes(x=x,y=y),
               fill=AreaColour2,colour=AreaLineColour)+
  geom_polygon(data=Cav3_25_2mgWT,mapping=aes(x=WT,y=Area90),
               fill=AreaColour1,colour=AreaLineColour)+
  geom_point(data=SumEBEsP3252[Age<=48],mapping=aes(x=WT,y=Cavss),
             colour=MarkerColour,size=3,shape=16,stroke=0.8)+
  geom_path(data=Cav3_25_2mgWT,mapping=aes(x=WT,y=Median),
            colour=LineColour,linewidth=1.2)+
  geom_vline(xintercept=10)+
  geom_vline(xintercept=20)+
  scale_x_continuous(breaks=seq(2,20,2))+
  coord_cartesian(ylim=c(0.5,5.3),xlim=c(2.9,21.4),expand=FALSE)+
  labs(y=expression('Drug X C'[av]*' (mg/L)'),
       x="Weight (kg)",
       title="Proposed posology")+
  themeO

pg4WT

print(ps2,position=c(0,0.5,0.5,1),more=TRUE)
print(ps2WT,position=c(0.5,0.5,1,1),more=TRUE)
print(ps4,position=c(0,0,0.5,0.5),more=TRUE)
print(ps4WT,position=c(0.5,0,1,0.5),more=FALSE)

(pg2+pg2WT)/(pg4+pg4WT)







