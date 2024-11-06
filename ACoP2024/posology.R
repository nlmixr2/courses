## posologyR

library(rxode2)
library(tidyverse)
library(monolix2rx)
library(posologyr)

# plot data

dat <- read_csv("monolix_data_acop/dapto_data_acop.csv")

ipreds <- read_csv("monolix_data_acop/daptomycin/IndividualParameters/estimatedIndividualParameters.txt")
ietas  <- read_csv("monolix_data_acop/daptomycin/IndividualParameters/estimatedRandomEffects.txt")

# Using the development version of monolix2rx you can produce the same
# model by:
# rx <- monolix2rx("monolix_data_acop/daptomycin.mlxtran")
#
# And then adding the covariates by piping
# rx <- rx %>% ini(SEX=1,WT=82.8, ClCr=115)

# model

rx <- function() {
  ini({
    Cl_pop <- -0.438992472277305
    Q_pop <- 1.29192666319929
    V1_pop <- 1.61789234074747
    V2_pop <- 1.63337359455679
    beta_Cl_SEX_1 <- -0.281709481904363
    beta_Cl_logtClCr <- 0.272203552838281
    beta_V2_logtWT <- 1.31586804217699
    SEX  <- 1
    WT   <- 82.8
    ClCr <- 115
    a <- 2.76801128307344
    omega_Cl ~ 0.0917131731665298
    omega_Q ~ 0.649552060001421
    omega_V1 ~ 0.231371471336246
    omega_V2 ~ 0.034723477760519
  })
  model({
    cmt(central)
    cmt(cmt2)

    Cl <- exp(Cl_pop + beta_Cl_SEX_1 * SEX + beta_Cl_logtClCr *
                log(ClCr/37.7806) + omega_Cl)
    Q <- exp(Q_pop + omega_Q)
    V1 <- exp(V1_pop + omega_V1)
    V2 <- exp(V2_pop + beta_V2_logtWT * log(WT/69.7821) +
                omega_V2)

    TVCL <- Cl
    TVV1 <- V1
    TVV2 <- V2
    TVQ <- Q
    k12 <- TVQ/TVV1
    k21 <- TVQ/TVV2

    d/dt(central) <- -k12 * central + k21 * cmt2 - TVCL/TVV1 *
      central

    d/dt(cmt2) <- +k12 * central - k21 * cmt2

    Cc <- central/TVV1
    DV <- Cc
    DV ~ add(a)
  })
}


# individual curves ############################################################

concs <- tibble()

getCurve <- function(id) {
  curdat    <- subset(dat, ID==id)
  curdos    <- subset(curdat, !is.na(AMT))
  curobs    <- subset(curdat, is.na(AMT))

  curid <- id
  curparams <- ipreds[ipreds$id==id,]

  parms <- tibble(
    Cl_pop = log(curparams$Cl_SAEM),
    V1_pop = log(curparams$V1_SAEM),
    Q_pop  = log(curparams$Q_SAEM),
    V2_pop = log(curparams$V2_SAEM),
    beta_Cl_SEX_1 <- 0,
    beta_Cl_logtClCr <- 0,
    beta_V2_logtWT <- 0,
    a <- 0,
    omega_Cl = 0,
    omega_Q = 0,
    omega_V1 = 0,
    omega_V2 = 0
  )

  et   <- et(id=id, amt=curdos$AMT, cmt=1, time=curdos$TIME, dur=curdos$DUR) %>%
    add.sampling(seq(0, max(curdat$TIME), length.out=500))

  crv <- rxSolve(rx, params = parms, events = et)
  crv$ID <- id

  concs <<- bind_rows(concs, crv)

}

for(i in unique(dat$ID)) {
  getCurve(i)
}

concs$time <- as.numeric(concs$time)
concs$Cc   <- as.numeric(concs$Cc)

ids <- sort(unique(concs$ID))

idr1 <- ids[1:8]
idr2 <- ids[9:16]
idr3 <- ids[17:24]
idr4 <- ids[25:32]
idr5 <- ids[33:40]
idr6 <- ids[41:48]
idr7 <- ids[49:56]
idr8 <- ids[57:64]
idr9 <- ids[65:72]

iplot <- function(idr) {
ggplot(concs[concs$ID %in% idr,], aes(as.numeric(time), as.numeric(Cc))) +
  geom_line() +
  geom_point(data=dat[dat$ID %in% idr & is.na(dat$AMT),], aes(TIME, DV), color="red") +
  scale_x_continuous("Time (h)") +
  scale_y_continuous("Concentration") +
  facet_wrap(~ factor(ID), nrow=2, scales = "free") +
  theme_light()

}

p1 <- iplot(idr1)
p2 <- iplot(idr2)
p3 <- iplot(idr3)
p4 <- iplot(idr4)
p5 <- iplot(idr5)
p6 <- iplot(idr6)
p7 <- iplot(idr7)
p8 <- iplot(idr8)
p9 <- iplot(idr9)

## examples ####################################################################

## patient 1

pat1 <- tibble(
  ID = 1,
  Cl_pop = -0.438992472277305,
  Q_pop = 1.29192666319929,
  V1_pop = 1.61789234074747,
  V2_pop = 1.63337359455679,
  beta_Cl_SEX_1 = -0.281709481904363,
  beta_Cl_logtClCr = 0.272203552838281,
  beta_V2_logtWT = 1.31586804217699,
  a = 0,
  omega_Cl = ietas$eta_Cl_SAEM[ipreds$id==1],
  omega_Q = ietas$eta_Q_SAEM[ipreds$id==1],
  omega_V1 = ietas$eta_V1_SAEM[ipreds$id==1],
  omega_V2 = ietas$eta_V2_SAEM[ipreds$id==1],
  SEX = ipreds$SEX[ipreds$id==1],
  WT = ipreds$WT[ipreds$id==1],
  ClCr = ipreds$ClCr[ipreds$id==1]
)

poso_dose_conc(dat = subset(dat, ID==1), prior_model = rx, time_c=6, time_dose=0, target_conc = 25, duration=0.25, indiv_param=pat1)

poso_dose_auc(dat = subset(dat, ID==1), prior_model = rx, starting_time = 168+144, time_auc=168+168, target_auc = 1000, duration=0.25, indiv_param=pat1, interdose_interval = 168, add_dose = 6)

poso_inter_cmin(dat = subset(dat, ID==1), prior_model = rx, duration=0.25, indiv_param=pat1, dose=1500, target_cmin=0.25)

poso_time_cmin(dat = subset(dat, ID==1), prior_model = rx, duration=0.25, indiv_param=pat1, dose=1500, target_cmin=25, from=0.5)

poso_estim_map(dat = subset(dat, ID==2), prior_model = rx)


## patient 17

pt17 <- ggplot(subset(dat, ID==17), aes(TIME, DV)) +
  geom_ribbon(data=tibble(TIME=c(-100,500), min=63.45, max=76.29, DV=0), aes(x=TIME, ymin=min, ymax=max), fill="blue", alpha=0.1) +
  geom_point(data=subset(dat, ID==17 & !is.na(DV)), col="red", size=3) +
  geom_vline(data=subset(dat, ID==17 & !is.na(AMT)), aes(xintercept=TIME), col="grey", linetype=2) +
  geom_text(data=subset(dat, ID==17 & !is.na(AMT)), aes(x=TIME, y=-5, label=paste(AMT, "mg")), angle=0, nudge_x = 0) +
  geom_hline(yintercept=c(63.45,76.29), col="blue", linetype=1) +
  geom_hline(yintercept=24.3, col="orange", linetype=1) +
  geom_hline(yintercept=0) +
  annotate("text", x=0, y=26.5, label="Cmin: toxicity cutoff of 24.3 mg/L", hjust=0, vjust=0, col="orange") +
  annotate("text", x=0, y=65.5, label="Cmax: efficacy range of 63.45-76.28 mg/L", hjust=0, vjust=0, col="blue") +
  scale_x_continuous("Time (h)", breaks=seq(0,500,by=24)) +
  scale_y_continuous("Concentration (mg/L)") +
  coord_cartesian(ylim=c(-10,175), xlim=c(0,192)) +
  theme_light() + theme(panel.grid = element_blank())

## now fit #####################################################################

pt17fit <- poso_estim_map(dat = subset(dat, ID==17), prior_model = rx)

pt17f <- pt17 +
  geom_line(data=pt17fit$model, aes(time, Cc), col="red")

## Cmax 70x the MIC ############################################################

pat17 <- tibble(
  ID = 17,
  Cl_pop = -0.438992472277305,
  Q_pop = 1.29192666319929,
  V1_pop = 1.61789234074747,
  V2_pop = 1.63337359455679,
  beta_Cl_SEX_1 = -0.281709481904363,
  beta_Cl_logtClCr = 0.272203552838281,
  beta_V2_logtWT = 1.31586804217699,
  a = 0,
  omega_Cl = pt17fit$model$params$omega_Cl,
  omega_Q = pt17fit$model$params$omega_Q,
  omega_V1 = pt17fit$model$params$omega_V1,
  omega_V2 = pt17fit$model$params$omega_V2,
  SEX = ipreds$SEX[ipreds$id==17],
  WT = ipreds$WT[ipreds$id==17],
  ClCr = ipreds$ClCr[ipreds$id==17]
)

## target Cmax of 70 mg/L
poso_dose_conc(dat = subset(dat, ID==17), prior_model = rx,
               time_c=(6*24)+0.5, time_dose=0,
               target_conc = 70,
               duration=0.5, indiv_param=pat17, add_dose=7,
               interdose_interval=24)

# $dose
# [1] 362.5421
#
# $type_of_estimate
# [1] "point estimate"
#
# $conc_estimate
# [1] 70
#
# $indiv_param
# # A tibble: 1 × 16
# ID Cl_pop Q_pop V1_pop V2_pop beta_Cl_SEX_1 beta_Cl_logtClCr beta_V2_logtWT     a omega_Cl
# <dbl>  <dbl> <dbl>  <dbl>  <dbl>         <dbl>            <dbl>          <dbl> <dbl>    <dbl>
#   1    17 -0.439  1.29   1.62   1.63        -0.282            0.272           1.32     0    0.258
# # ℹ 6 more variables: omega_Q <dbl>, omega_V1 <dbl>, omega_V2 <dbl>, SEX <dbl>, WT <dbl>,
# #   ClCr <dbl>

## use the dose we just calculated

et <- et() %>%
  add.dosing(dose=362.5421, nbr.doses = 7, dosing.interval = 24, rate = 725.069) %>%
  add.sampling(seq(0, 192, length.out=2000)) %>%
  add.sampling(c(0.5,24.5,48.5,72.5,96.5,120.5,144.5,168.5))

rx17 <- rxSolve(rx, params = pat17, events = et)

pt17optdose <- ggplot(subset(dat, ID==17), aes(TIME, DV)) +
  geom_ribbon(data=tibble(TIME=c(-100,500), min=63.45, max=76.29, DV=0), aes(x=TIME, ymin=min, ymax=max), fill="blue", alpha=0.1) +
  geom_vline(xintercept=c(0,24,48,72,96,120,144), col="grey", linetype=2) +
  geom_text(data=tibble(TIME=c(0,24,48,72,96,120,144)), aes(x=TIME, y=-5, label="362.5 mg"), angle=0, nudge_x = 0) +
  geom_hline(yintercept=c(63.45,76.29), col="blue", linetype=1) +
  geom_hline(yintercept=70, col="blue", linetype=2) +
  geom_hline(yintercept=24.3, col="orange", linetype=1) +
  geom_hline(yintercept=0) +
  geom_line(data=rx17, aes(as.numeric(time), as.numeric(Cc)), col="red") +
  annotate("text", x=0, y=26.5, label="Cmin: toxicity cutoff of 24.3 mg/L", hjust=0, vjust=0, col="orange") +
  annotate("text", x=0, y=65.5, label="Cmax: efficacy range of 63.45-76.28 mg/L", hjust=0, vjust=0, col="blue") +
  annotate("text", x=0, y=71.5, label="Cmax target: 70 mg/L", hjust=0, vjust=0, col="blue") +
  scale_x_continuous("Time (h)", breaks=seq(0,500,by=24)) +
  scale_y_continuous("Concentration (mg/L)") +
  coord_cartesian(ylim=c(-10,80), xlim=c(0,168)) +
  theme_light() + theme(panel.grid = element_blank())


## AUC of 700 ##################################################################

poso_dose_auc(dat = subset(dat, ID==17),
              prior_model = rx, starting_time = 120,
              time_auc=24,
              target_auc = 700,
              duration=0.5, indiv_param=pat17,
              interdose_interval = 24, add_dose = 7)

# $dose
# [1] 597.3413
#
# $type_of_estimate
# [1] "point estimate"
#
# $auc_estimate
# [1] 700
#
# $indiv_param
# # A tibble: 1 × 16
# ID Cl_pop Q_pop V1_pop V2_pop beta_Cl_SEX_1 beta_Cl_logtClCr beta_V2_logtWT     a omega_Cl
# <dbl>  <dbl> <dbl>  <dbl>  <dbl>         <dbl>            <dbl>          <dbl> <dbl>    <dbl>
#   1    17 -0.439  1.29   1.62   1.63        -0.282            0.272           1.32     0    0.258
# # ℹ 6 more variables: omega_Q <dbl>, omega_V1 <dbl>, omega_V2 <dbl>, SEX <dbl>, WT <dbl>,
# #   ClCr <dbl>

## use the dose we just calculated

et <- et() %>%
  add.dosing(dose=597.3413, nbr.doses = 7, dosing.interval = 24, rate = 1194.647) %>%
  add.sampling(seq(0, 192, length.out=2000)) %>%
  add.sampling(c(0.5,24.5,48.5,72.5,96.5,120.5,144.5,168.5))

rx17 <- rxSolve(rx, params = pat17, events = et)


pol <- tibble(
  x = c(rx17$time[rx17$time>=120 & rx17$time <144], rev(rx17$time[rx17$time>=120 & rx17$time <144])),
  y = c(rx17$Cc[rx17$time>=120 & rx17$time <144], rep(0, times=length(rx17$time[rx17$time>=120 & rx17$time <144]))))


pt17optauc <- ggplot(subset(dat, ID==17), aes(TIME, DV)) +
  geom_ribbon(data=tibble(TIME=c(-100,500), min=63.45, max=76.29, DV=0), aes(x=TIME, ymin=min, ymax=max), fill="blue", alpha=0.1) +
  geom_vline(xintercept=c(0,24,48,72,96,120,144), col="grey", linetype=2) +
  geom_text(data=tibble(TIME=c(0,24,48,72,96,120,144)), aes(x=TIME, y=-5, label="597.3 mg"), angle=0, nudge_x = 0) +
  geom_hline(yintercept=c(63.45,76.29), col="blue", linetype=1) +
  geom_hline(yintercept=70, col="blue", linetype=2) +
  geom_hline(yintercept=24.3, col="orange", linetype=1) +
  geom_hline(yintercept=0) +
  geom_line(data=rx17, aes(as.numeric(time), as.numeric(Cc)), col="red") +
  geom_polygon(data=pol, aes(x=x, y=y), fill="red", alpha=0.1) +
  annotate("text", x=0, y=26.5, label="Cmin: toxicity cutoff of 24.3 mg/L", hjust=0, vjust=0, col="orange") +
  annotate("text", x=0, y=65.5, label="Cmax: efficacy range of 63.45-76.28 mg/L", hjust=0, vjust=0, col="blue") +
  annotate("text", x=122, y=5, label="AUC: 700 mg.d/L", hjust=0, vjust=0, col="red") +
  scale_x_continuous("Time (h)", breaks=seq(0,500,by=24)) +
  scale_y_continuous("Concentration (mg/L)") +
  coord_cartesian(ylim=c(-10,125), xlim=c(0,168)) +
  theme_light() + theme(panel.grid = element_blank())


## simulate using this subject's EBEs ##########################################

sdat <- data.frame(ID=17,
                   TIME=seq(0, 192, length.out=2000),
                   DV=0,
                   AMT=0,
                   EVID=0,
                   DUR=NA)

sdos <- data.frame(ID=17,
                   TIME=seq(0, 144, by=24),
                   DV=NA,
                   AMT=750,
                   EVID=1,
                   DUR=0.5)

sdats <- rbind(sdat, sdos)
sdats <- sdats[order(sdats$TIME),]

# run the simulations

sim17 <- poso_simu_pop(dat=sdats,prior_model=rx,n_simul=2000, return_model=T)

sim17s <- sim17$model %>%
  dplyr::group_by(time) %>%
  dplyr::summarise(mean=median(Cc),
                   lower50=quantile(Cc,0.25),
                   upper50=quantile(Cc,0.75),
                   lower90=quantile(Cc,0.05),
                   upper90=quantile(Cc,0.95),
                   lower95=quantile(Cc,0.025),
                   upper95=quantile(Cc,0.975))

sim17s$DV <- 0

# plot them

plotsim17 <- ggplot(subset(dat, ID==17), aes(TIME, DV)) +
  geom_ribbon(data=tibble(TIME=c(-100,500), min=63.45, max=76.29, DV=0), aes(x=TIME, ymin=min, ymax=max), fill="blue", alpha=0.1) +
  geom_ribbon(data=sim17s, aes(x=time, ymin=lower95, ymax=upper95), fill="red", alpha=0.1) +
  geom_ribbon(data=sim17s, aes(x=time, ymin=lower90, ymax=upper90), fill="red", alpha=0.1) +
  geom_ribbon(data=sim17s, aes(x=time, ymin=lower50, ymax=upper50), fill="red", alpha=0.1) +
  geom_line(data=sim17s, aes(x=as.numeric(time), y=mean), col="red") +
  geom_text(data=tibble(TIME=c(0,24,48,72,96,120,144)), aes(x=TIME, y=-5, label="750 mg"), angle=0, nudge_x = 0) +
  geom_vline(xintercept=c(0,24,48,72,96,120,144), col="grey", linetype=2) +
  geom_hline(yintercept=c(63.45,76.29), col="blue", linetype=1) +
  geom_hline(yintercept=24.3, col="orange", linetype=1) +
  geom_hline(yintercept=0) +
  annotate("text", x=5, y=28.5, label="Cmin: toxicity cutoff of 24.3 mg/L", hjust=0, vjust=0, col="black") +
  annotate("text", x=5, y=67.5, label="Cmax: efficacy range of 63.45-76.28 mg/L", hjust=0, vjust=0, col="blue") +
  scale_x_continuous("Time (h)", breaks=seq(0,500,by=24)) +
  scale_y_continuous("Concentration (mg/L)") +
  coord_cartesian(ylim=c(-10,250), xlim=c(0,168))+
  theme_light() + theme(panel.grid = element_blank())
