## load the required libraries
library(nlmixr2)
library(xpose.nlmixr2)
library(ggPMX)
library(xgxr)
library(tidyverse)
library(ggplot2)
library(patchwork)




#################################################################################
##                                                                             ##
## nlmixr analysis Part 1                                                      ##
##                                                                             ##
## Warfarin population PK using SAEM estimation                                ##
##                                                                             ##
#################################################################################

## Use the Warfarin data from nlmixr and subset to PK-only using tidyverse
PKdata <- nlmixr2data::warfarin %>%
  filter(dvid=="cp") %>% # use dplyr to filter to "cp" concentrations
  select(-dvid) %>% # remove dvid from the dataset
  mutate(sparse=factor(ifelse(time < 24, "Day 1 (Intense)", "After Day 1 (Sparse)"), # create a flag for sparse PK
                       levels=c("Day 1 (Intense)", "After Day 1 (Sparse)")))  # Set the order for intense followed by Sparse
  
## look at data structure and identify what each column is for
#View(PKdata)


## Our warfarin data file: a simple ggplot2 plot to provide an impression


# as with all ggplots you can comment out a line to see what it does...
gg <- ggplot(PKdata, aes(time, dv, color=factor(id))) + # declare a dv vs time plot grouped by ID
  geom_line(size=1.1) +   # Add lines with a slightly thicker color
  geom_point() +  # Add points to the plot
  xgx_theme() + # Make the plot follow the graphics principles
                # (https://opensource.nibr.com/xgx/Resources/Graphics_Principles_Cheat_Sheet_v1.1.pdf)
  xgx_scale_x_time_units(units_dataset = "hr", units_plot = "hr") + # display tick marks in more time-friendly manner
  xlab("Time (h)") + # X axis label
  ylab("Warfarin (mg/L)") + # Y axis label
  theme(legend.position = "none") # remove legend

print(gg)

# What does this mean? What if we stratify by the sparse vs the full pk?

gg2 <- gg + facet_wrap(~sparse, scales="free_x")

print(gg2)

# We can also see this on a semi-log plot
gg3 <- gg2 +
  xgx_scale_y_log10()

print(gg3)

# What does this mean?

# We briefly used some exploratory analysis using ggplot and xgx helper functions

## Use xgxr for simplified concentration over time, colored by Dose, mean +/- 95% CI

#Often in exploring data it is worthwhile to plot by dose by each nominal
#time and add the 95% confidence interval. This typical plot can be
#cumbersome and lack some nice features that `xgxr` can help with. Note
#the following helper functions:
#
#-   `xgx_theme_set()` this sets the theme to black and white color theme
#    and other best pratices in `xgxr`.
#
#-   `xgx_geom_ci()` which creates the Confidence Interval and mean plots
#    in a simple interface.
#
#-   `xgx_scale_y_log10()` which creates a log-scale that includes the
#    minor grids that immediately show the viewer that the plot is a
#    semi-log plot without carefully examining the y axis.
#
#-   `xgx_scale_x_time_units()` which creates an appropriate scale based
#    on your times observed and the units you use. It also allows you to
#    convert units easily for the right display.
#
#-   `xgx_annote_status()` which adds a `DRAFT` annotation which is often
#    considered best practice when the data or plots are draft.


#Not only is it useful to look at the mean concentrations, it is often
#useful to look at the mean concentrations and their relationship between
#actual individual profiles. Using `ggplot` coupled with the `xgxr`
#helper functions used above, we can easily create these plots as well:

xgx_theme_set()

## Exploratory xgx graph on logarithmic scale with mean and SD error bars at scheduled time points
gg1 <- ggplot(data = PKdata, aes(x = time, y = dv)) +
  geom_line(aes(group = id), color = "grey50", size = 1, alpha = 0.3) +
  xgx_geom_ci(aes(x = time, color = NULL, group = NULL, shape = NULL), conf_level = 0.95) +
  xgx_scale_y_log10() +
  xgx_scale_x_time_units(units_dataset = "hours", units_plot = "days") +
  labs(y = "Warfarin (mg/L)") +
  theme(legend.position = "none")

print(gg1)

# As an exercise you could focus on the first part of the data or with change the type of confidence interval 


# Next, define a first order-absorption linear elimination model using a solved system solution

One.comp.KA.solved <- function() {
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
    ## solved system example
    ## where residual error is assumed to follow proportional and additive error
    linCmt() ~ prop(prop.err) + add(add.err)
  })
}


## clean the cache to be sure there are no remaining old models
rxClean()

## Check the model and some of the assumptions made by nlmixr 
nlmixr(One.comp.KA.solved) 

## estimate parameters using nlmixr:
fitOne.comp.KA.solved_S <-
  nlmixr(
    One.comp.KA.solved,    #the model definition
    PKdata,                #the data set
    est = "saem",          #the estimation algorithm (SAEM)
                             #the SAEM minimisation options:
    saemControl(nBurn = 200, #200 SAEM burn-in iterations (the default)
                nEm   = 300, #300 EM iterations (the default)
                print = 50), #print every 50th iteration
    #only print every 50th estimation step (default=1 which gives endless output)
    tableControl(cwres = TRUE,npde=TRUE) #calculates NONMEM-style conditional weighted residuals and npde for diagnostics
  )

## results are stored in the nlmixr object and can be viewed:
fitOne.comp.KA.solved_S

## for SAEM, convergence can be checked using a parameter trace plot:
traceplot(fitOne.comp.KA.solved_S)

# Model Diagnostics with ggPMX
# The controller is first constructed
ctr <- pmx_nlmixr(fitOne.comp.KA.solved_S, 
                  conts = c("wt","age"),
                  cats=c("sex"), 
                  vpc=FALSE,
                  settings=pmx_settings(is.draft=FALSE))

# Printing the object lists the plots that can be generated
ctr

# and can then be piped into a specific plot
# syntax for npde vs pred plot
ctr %>% pmx_plot_npde_pred

# alternatively:
pmx_plot_npde_pred(ctr)

## Modify graphical options and add DRAFT label:
ctr %>% pmx_plot_npde_time(smooth = list(color="blue"), 
                           point = list(shape=4), 
                           s.draft=TRUE, 
                           labels = list(x = "Time (hours)", 
                                         y = "Normalized PDE"))

## DV vs IPRED plot
ctr %>% pmx_plot_dv_ipred(scale_x_log10=TRUE, scale_y_log10=TRUE)

#You can filter to restrict the values of IPRED for instance:
ctr %>% pmx_plot_dv_ipred(scale_x_log10=TRUE, scale_y_log10=TRUE,filter=(IPRED>1))

## DV vs PRED plot
ctr %>% pmx_plot_dv_pred(scale_x_log10=TRUE, scale_y_log10=TRUE)

## Absolute individual weighted residuals to investigate the residual error model
ctr %>% pmx_plot_abs_iwres_ipred
#again, alternatively:
#pmx_plot_abs_iwres_ipred(ctr)

## Absolute individual weighted residuals vs time to detect model discrepancies
pmx_plot_abs_iwres_time(ctr)

ctr %>% pmx_plot_iwres_dens
ctr %>% pmx_plot_eta_qq
ctr %>% pmx_plot_eta_box
ctr %>% pmx_plot_eta_hist
ctr %>% pmx_plot_eta_matrix

## generate a full report with diagnostics

ctr %>% pmx_report(name="ggPMX_report",
                   save_dir=".",
                   format="report",
                   extension="word")


## Instead of using ggPMX, an alternative is to use the new xpose package for which xpose.nlmixr provides the link
## the nlmixr object can be transformed into an xpose object to allow diagnostics with the new xpose package
## the link between nlmixr and xpose is provided by the xpose.nlmixr package
## only xpose_data_nlmixr is from xpose.nlmixr, all further commands are from the xpose package
xpdb.1s <- xpose_data_nlmixr(fitOne.comp.KA.solved_S)

## this can also be used to generate trace plots (parameters vs iterations:)
prm_vs_iteration(xpdb.1s)
## to remove the path to the script from the plot use:
prm_vs_iteration(xpdb.1s,caption=NULL)

## and many types of diagnostic plots (see cheatsheet)

## dv vs cpred plot:
dv_vs_pred(xpdb.1s,
           caption = NULL)

# by default model typical predictions (PRED) are assigned to CPRED (conditional population predictions):
list_vars(xpdb.1s)
# if you want this to be PRED instead, these can be updated using 'standard' syntax:
xpdb.1s<-set_var_types(xpdb.1s,pred = 'PRED')
# or using magrittr piping type code:
# xpdb.1s<-xpdb.1s %>% set_var_types(pred = 'PRED')

## dv vs pred plot:
dv_vs_pred(xpdb.1s,
           caption = NULL)

## dv vs ipred plot:
dv_vs_ipred(xpdb.1s,        #the xpose object
            caption = NULL) #if not NULL provides the directory where this was run

## CWRES vs time:
res_vs_idv(xpdb.1s,           #the xpose object
           res = "CWRES",     #examine conditional weighted residuals
           idv = "TIME",      #as a function of time
           caption = NULL)    #if not NULL provides the directory where this was run

#Absolute values of individual weighted residual vs time
absval_res_vs_idv(xpdb.1s,        #the xpose object
                  res = "IWRES",  #examine absolute values (absval) of individual weighted residuals
                  idv = "TIME",   #as a function of time
                  caption = NULL) #if not NULL provides the directory where this was run

## ...this plot shows that there are some issues at the earlier time points that an updated absorption model may fix



#######################################################################
###                                                                 ###
###   Visual predictive checks                                      ###
###                                                                 ###
#######################################################################


## nlmixr comes with its built-in vpc functionality that uses Ron Keizer's vpc package
## see the cheatsheet for further options

## because the data set uses nominal time points, it is nice to have the bins surround these time points 
## so that each time point falls in a bin

bin_mids <- sort(unique(PKdata$time))
bin_edges <- bin_mids - c(0, diff(bin_mids) / 2) # puts the edges in the middle of the nominal time points

p1<-vpcPlot(
  fitOne.comp.KA.solved_S,        #the nlmixr object
  n = 500,                        #number of simulated trials
  bins = bin_edges,
  show = list(                    #additional items to show, like the observations
    obs_dv = TRUE,
    obs_median = TRUE,
    sim_median = TRUE,
    sim_median_ci = TRUE,
    obs_ci = TRUE,
    pi = TRUE
  ), 
  xlab = "Time (h)",              #x-axis label
  ylab = "Concentration (mg/L)",  #y-axis label
  title = "VPC for first order absorption PopPK model\nwith linear y axis"
)

print(p1)



## or with a log y-axis starting at 0.5
p2<-vpcPlot(
  fitOne.comp.KA.solved_S,
  n = 500,
  bins = bin_edges,
  show = list(                    #additional items to show, like the observations
    obs_dv = TRUE,
    obs_median = TRUE,
    sim_median = TRUE,
    sim_median_ci = TRUE,
    obs_ci = TRUE,
    pi = TRUE), 
  xlab = "Time (h)",
  ylab = "Concentration (mg/L)",
  title = "VPC for first order absorption PopPK model\nwith log y-axis",
  log_y = TRUE,            #to request a log y-axis
  log_y_min = 0.5          #starting at 0.5
)
print(p2)

#Using patchwork syntax:
print(p1+p2)



#Alternatively using direct Ron Keizer code so you only need to simulate once:
library(vpc)

SimVPC<-vpcSim(
  fitOne.comp.KA.solved_S,                  #the nlmixr object
  n = 500)


# rename the sim to dv

SimVPC <- SimVPC %>%
  rename(dv=sim)


p0 <- vpc_vpc(sim = SimVPC, obs = PKdata,bins=bin_edges,
        show = list(obs_dv = TRUE,
                    obs_median=TRUE,
                    sim_median=TRUE,
                    sim_median_ci=TRUE,
                    obs_ci=TRUE,
                    pi=TRUE),
        xlab = "Time (h)",              #x-axis label
        ylab = "Concentration (mg/L)",  #y-axis label
        title = "VPC for first order absorption PopPK model\nwith linear y axis")

print(p0)



# ## Individual fits can be generated using using xpose:
ind_plots(xpdb.1s,caption = NULL,ncol = 4,nrow = 4)

# ## ...use the arrows in the plot window to examine the earlier curves

## Individual fits can also be generated using augPred (augmented predictions)
## that provides smooth profiles by interpolating the predictions between observations:
plot(augPred(fitOne.comp.KA.solved_S))
## ...use the arrows in the plot window to examine the earlier curves

#or the augPred output can be formatted to your liking for instance using whatever tool you like, since it is data
indivpk<-augPred(fitOne.comp.KA.solved_S)

# Exercise: create individual plots the way you wish using the augPred dataset.

#################################################################################
##                                                                             ##
## nlmixr analysis Part 2                                                      ##
##                                                                             ##
## Warfarin population PK using FOCEi and ODEs                                 ##
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

rxClean()
## estimate parameters using nlmixr/FOCEI:
fitOne.comp.KA.ODE_F <-
  nlmixr(One.comp.KA.ODE,          #the model definition
         PKdata,                   #the data set
         est = "focei",            #the estimation algorithm (FOCEi)
                                   #FOCEi options:
         foceiControl(print = 5))  #only print every 5th estimation step

ctr <- pmx_nlmixr(fitOne.comp.KA.ODE_F, conts = c("wt","age"),cats=c("sex"),
                  vpc=FALSE,settings=pmx_settings(is.draft=FALSE))



ctr %>% pmx_plot_npde_pred
ctr %>% pmx_plot_npd_pred
ctr %>% pmx_plot_cwres_pred

ctr %>% pmx_report(name="ggPMX_report2",
                   save_dir=".",
                   format="report",
                   extension="word")

#################################################################################
##                                                                             ##
##  Hands-on assignment: nlmixr model development                              ##
##  Examine the GOF plots and implement alternative absorption models          ##
##  -one or more transit compartment(s)                                        ##
##                                                                             ##
##  Compare vpcs of alternatives and compare OFVs:                             ##
##                                                                             ##
##  fitPK001$OBJF-fitPK002$OBJF                                                ##
##                                                                             ##
#################################################################################

#...

#################################################################################
##                                                                             ##
#################################################################################

