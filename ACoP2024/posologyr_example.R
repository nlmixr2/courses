
# AUC dose selection example

library(posologyr)

# AUC-based dosage adjustment for a patient treated with vancomycin
# for methicillin-resistant Staphylococcus aureus blood stream
# infection, using the population pharmacokinetic (ppk) model of Goti
# et al. 2018, using the data from therapeutic drug monitoring (TDM).


mod_vancomycin_Goti2018 <- function() {
  ini({
    THETA_Cl <- 4.5
    THETA_Vc <- 58.4
    THETA_Vp <- 38.4
    THETA_Q <- 6.5
    ETA_Cl ~ 0.147
    ETA_Vc ~ 0.510
    ETA_Vp ~ 0.282
    add.sd <- 3.4
    prop.sd <- 0.227
  })
  model({
    TVCl  = THETA_Cl*(CLCREAT/120)^0.8*(0.7^DIAL);
    TVVc  = THETA_Vc*(WT/70)          *(0.5^DIAL);
    TVVp  = THETA_Vp;
    TVQ   = THETA_Q;
    Cl    = TVCl*exp(ETA_Cl);
    Vc    = TVVc*exp(ETA_Vc);
    Vp    = TVVp*exp(ETA_Vp);
    Q     = TVQ;
    ke    = Cl/Vc;
    k12   = Q/Vc;
    k21   = Q/Vp;
    Cc    = centr/Vc;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC) <- Cc # NOTE AUC is calculated in the model
    Cc ~ add(add.sd) + prop(prop.sd) + combined1()
  })
}


# Discontinuous intravenous infusion

df_patientB <- data.frame(ID=1,TIME=c(0.0,13.0,24.2,48),
                          DV=c(NA,12,NA,9.5),
                          AMT=c(2000,0,1000,0),
                          DUR=c(2,NA,2,NA),
                          EVID=c(1,0,1,0),
                          CLCREAT=65,WT=70,DIAL=0)
df_patientB

# Estimate the MAP individual parameters

patB_map <- poso_estim_map(dat=df_patientB,
                           prior_model=mod_vancomycin_Goti2018)

#Plot the individual pharmacokinetic profile

# The individual pharmacokinetic profile can be plotted using the
# rxode2 model provided by the poso_estim_map() function.

plot(patB_map$model,Cc)

# Using ggplot2 the observed data points can be added to the plot

library(ggplot2)

#Get the observations from the patient record
indiv_obs             <- df_patientB[,c("DV","TIME")]
names(indiv_obs)      <- c("value","time")

#Overlay the MAP profile and the observations
plot(patB_map$model,Cc) +
  ylab("Central concentration") +
  geom_point(data=indiv_obs, size= 3, na.rm=TRUE)

# The MAP profile matches the observations.

# Get the AUC24 from the MAP model

# Considering a MIC of 1 mg/L, the target AUC over 24 hours (AUC24) is
# 400 mg.h/L. The AUC can be retrieved from the rxode2 model using the
# usual R data.frame syntax.

#AUC 0_24
AUC_map_first_dose <- patB_map$model$AUC[which(patB_map$model$time == 24)]
AUC_map_first_dose

#AUC 24_48
AUC_map_second_dose <- patB_map$model$AUC[which(patB_map$model$time == 48)] - AUC_map_first_dose

AUC_map_second_dose

# The current dosage does not meet the target AUC.

##########################################################################
# Optimal maintenance dose selection a posteriori

# The maintenance dose needed to reliably achieve an AUC24 of 400
# mg*h/L can be estimated by simulating a multiple dose regimen over
# enough administrations (e.g. 11 consecutive administrations, with
# add_dose=10) to approximate the steady-state.

poso_dose_auc(dat=df_patientB,
              prior_model=mod_vancomycin_Goti2018,
              time_auc=24,
              starting_time=24*9,
              interdose_interval=24,
              add_dose=10,
              duration=2,
              target_auc=400)

# This gives an optimal dose of ~1200


# Continuous intravenous infusion

# The maintenance dose for a continuous intravenous infusion can be
# easily determined by setting the duration of the infusion equal to
# the interdose_interval.

poso_dose_auc(dat=df_patientB,
              prior_model=mod_vancomycin_Goti2018,
              time_auc=24,
              starting_time=24*9,
              interdose_interval=24,
              add_dose=10,
              duration=24,
              target_auc=400)


# The optimal maintenance dose is also 1200 mg / 24 h for a continuous intravenous infusion.
