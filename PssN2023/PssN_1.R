## load the required libraries
library(data.table)
library(xgxr)
library(ggplot2)
library(nlmixr2)
library(xpose.nlmixr2)
library(ggPMX)
library(patchwork)
library(lattice)



#################################################################################
##                                                                             ##
## nlmixr2 analysis Part 1                                                      ##
##                                                                             ##
## Warfarin population PK using SAEM estimation                                ##
##                                                                             ##
#################################################################################

## read in the Warfarin PK-only data set using data.table syntax (fast and efficient!)
PKdata <- fread("warfarin_PKS.csv")

## look at data structure and identify what each column is for
#View(PKdata)
PKdata[,Group:=factor(SPARSE,levels=c(0,1),labels=c("With day 1","From day 1 only"))]

PKdata[,sparse:=factor(ifelse(TIME < 24, "Day 1 (Intense)", "After Day 1 (Sparse)"), # create a flag for sparse PK
                       levels=c("Day 1 (Intense)", "After Day 1 (Sparse)"))]  # Set the order for intense followed by Sparse



#################################################################################
##                                                                             ##
## Exploratory analysis using ggplot and xgx helper functions                  ##
##                                                                             ##
#################################################################################

# Use xgxr for summarised concentration over time, split by time period
#
# Often in exploring data it is worthwhile to plot by a factor (like dose) 
# by each nominal time and add the 95% confidence interval. This typical plot
# can be cumbersome to create and lack some nice features that `xgxr` can help 
# with. Note the following helper functions:
#
# -   `xgx_theme_set()` this sets the theme to black and white color theme
#    and other best practices in `xgxr`.
#
# -   `xgx_geom_ci()` which creates the Confidence Interval and mean plots
#    in a simple interface.
#
# -   `xgx_scale_y_log10()` which creates a log-scale that includes the
#    minor grids that immediately show the viewer that the plot is a
#    semi-log plot without carefully examining the y axis.
#
# -   `xgx_scale_x_time_units()` which creates an appropriate scale based
#    on your times observed and the units you use. It also allows you to
#    convert units easily for the right display.
#
# -   `xgx_annote_status()` which adds a `DRAFT` annotation which is often
#    considered best practice when the data or plots are draft.


xgx_theme_set() # This uses black and white theme based on xgxr best practices


# Not only is it useful to look at the mean concentrations, it is often
# useful to look at the mean concentrations and their relationship between
# actual individual profiles. Using `ggplot` coupled with the `xgxr`
# helper functions used above, we can easily create these plots as well: 


ggp1<-ggplot(data = PKdata, aes(x = TIME, y = DV)) +
  geom_line(aes(group = ID), color = "grey50",
              size = 1, alpha = 0.3) +
  labs(y = "Warfarin (mg/L)") +
  theme(legend.position = "none") 
print(ggp1)

#Change the x-axis from hours to days and add a proper label

ggp2 <- ggp1 + xgx_scale_x_time_units(units_dataset = "hours", units_plot = "days") 
print(ggp2)

#If we stratify by the two groups in the data:

ggp3 <- ggp2 + facet_grid(.~Group) 
print(ggp3)

#Or alternatively split between rich Day 1 and later sparse data:

ggp4 <- ggp3 + facet_wrap(~sparse, scales="free_x")
print(ggp4)

#Because there are only two groups with information, the same is given using facet_grid in this case
ggp4 <- ggp3 + facet_grid(~sparse, scales="free_x")
print(ggp4)

# We can also see this on a semi-log plot
ggp5 <- ggp4 + xgx_scale_y_log10()
print(ggp5)

#Clues towards what model to use?

# xgx can also add nice summary information if data has similar sampling times (i.e., when nominal times are provided) 
# summaries of mean plus 95% CI
ggp6 <- ggp5 + xgx_geom_ci(aes(x = TIME, color = NULL, group = NULL, shape = NULL), conf_level = 0.95)
print(ggp6)

#On linear scale this would result in a CI crossing zero because CIs are assumed symmetrical:
ggp7 <- ggp4 + xgx_geom_ci(aes(x = TIME, color = NULL, group = NULL, shape = NULL), conf_level = 0.95)
print(ggp7)

#Perhaps a median and 95% of the data would be more suitable:
ggp8 <- ggp4 + xgx_geom_pi(aes(x = TIME, color = NULL, group = NULL, shape = NULL)) 
print(ggp8)



## Define a first order-absorption linear elimination model using a solved system solution

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

## Check the model and some of the assumptions made by nlmixr2 
nlmixr2(One.comp.KA.solved) 

## estimate parameters using nlmixr2:
run001S <-
  nlmixr2(
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

## results are stored in the nlmixr2 object and can be viewed:
run001S

## and saved for future use or reference:
save(run001S, file = "run001S.Rdata")
#load(file = "run001S.Rdata")

## for SAEM, convergence can be checked using a parameter trace plot:
traceplot(run001S)



# Model Diagnostics with ggPMX
# The controller is first constructed

ctr001S <- pmx_nlmixr(run001S, 
                  conts = c("WT","AGE"),
                  cats=c("SEX","SPARSE"), 
                  vpc=FALSE,
                  settings=pmx_settings(is.draft=FALSE))

# and can then be piped into a specific plot
# syntax for npde vs pred plot
ctr001S %>% pmx_plot_npde_pred
# alternatively:
pmx_plot_npde_pred(ctr001S)

## Alternative graph (NPDE) with modified graphical options
ctr001S %>% pmx_plot_npde_time(smooth = list(color="blue"), 
                           point = list(shape=4), 
                           is.draft=TRUE, 
                           labels = list(x = "Time (hours)", 
                                         y = "Normalized PDE"))

## DV vs IPRED plot
ctr001S %>% pmx_plot_dv_ipred(scale_x_log10=TRUE, scale_y_log10=TRUE)

#You can filter to restrict the values of IPRED for instance:
ctr001S %>% pmx_plot_dv_ipred(scale_x_log10=TRUE, scale_y_log10=TRUE,filter=(IPRED>1))

## DV vs PRED plot
ctr001S %>% pmx_plot_dv_pred(scale_x_log10=TRUE, scale_y_log10=TRUE)

## Absolute individual weighted residuals to investigate the residual error model
ctr001S %>% pmx_plot_abs_iwres_ipred

## Absolute individual weighted residuals vs time to detect model discrepancies
ctr001S %>% pmx_plot_abs_iwres_time

ctr001S %>% pmx_plot_iwres_dens
ctr001S %>% pmx_plot_eta_qq
ctr001S %>% pmx_plot_eta_box
ctr001S %>% pmx_plot_eta_hist
ctr001S %>% pmx_plot_eta_matrix

## generate a full report with diagnostics

ctr001S %>% pmx_report(name="ggPMX_report001S",
                   save_dir=".",
                   format="report",
                   extension="word")


## Instead of using ggPMX, an alternative is to use the new xpose package for which xpose.nlmixr2 provides the link
## the nlmix2 object can be transformed into an xpose object to allow diagnostics with the new xpose package
## the link between nlmixr2 and xpose is provided by the xpose.nlmixr2 package
## only xpose_data_nlmixr is from xpose.nlmixr, all further commands are from the xpose package
xpdb.1s <- xpose_data_nlmixr(run001S)

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


## nlmixr2 comes with its built-in vpc functionality that uses Ron Keizer's vpc package
## see the cheatsheet for further options

## because the data set uses nominal time points, it is nice to have the bins surround these time points 
## so that each time point falls in a bin

bin_mids <- sort(unique(PKdata$TIME))
bin_edges <- bin_mids - c(0, diff(bin_mids) / 2) # puts the edges in the middle of the nominal time points

p1<-vpcPlot(
  run001S,        #the nlmixr object
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
p1



## or with a log y-axis starting at 0.5
p2<-vpcPlot(
  run001S,
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
p2

#Using patchwork syntax:
p1+p2



#Alternatively using direct Ron Keizer code so you only need to simulate once:
library(vpc)

SimVPC<-vpcSim(
  run001S,                  #the nlmixr object
  n = 500)


setnames(SimVPC,c("id","time","sim"),c("ID","TIME","DV"))

vpc_vpc(sim = SimVPC, obs = PKdata,bins=bin_edges,
        show = list(obs_dv = TRUE,
                    obs_median=TRUE,
                    sim_median=TRUE,
                    sim_median_ci=TRUE,
                    obs_ci=TRUE,
                    pi=TRUE),
        xlab = "Time (h)",              #x-axis label
        ylab = "Concentration (mg/L)",  #y-axis label
        title = "VPC for first order absorption PopPK model\nwith linear y axis")





# ## Individual fits can be generated using using xpose:
ind_plots(xpdb.1s,caption = NULL,ncol = 4,nrow = 4)
# ## ...use the arrows in the plot window to examine the earlier curves

## Individual fits can also be generated using augPred (augmented predictions)
## that provides smooth profiles by interpolating the predictions between observations:
plot(augPred(run001S))
## ...use the arrows in the plot window to examine the earlier curves

#or the augPred output can be formatted to your liking for instance using lattice:
indivpk<-augPred(run001S)

table(indivpk$ind)
#Population Individual   Observed 
#      1883       1883        251

## specify array of colours for curves
nlmixCOLS <- c("#28466A","#8DB6CD","#B40000")  

xyplot(
  values~time|id,          ## plot the variable 'values' by time and make a separate panel for each id
  data=indivpk,            ## data source
  groups=ind,              ## make separate curves for ind that indicates Observed data, 
                           ## Indivdual predictions and Population predictions
  layout=c(8,4),           ## arrange as 8 columns and 4 rows
  type=c("l","l","p"),     ## represent these three by a line, a line and only markers (l=line, p=points)
  col=nlmixCOLS[c(2,1,3)], ## colours for each curve
  cex=c(0.1,0.1,1),        ## character size for the markers
  lwd=c(2,2,0.1),          ## line width of the lines
  pch=19,                  ## use closed circles as marker
  xlab="Time (hr)\n",      ## x-axis label
  ylab="Warfarin (mg/L)",  ## y-axis label
  as.table=TRUE,           ## have the first plot at the top left (otherwise plots start at the lower left corner)
  scales=list(alternating=1),  ## have axis labels at left and bottom (and not alternating)
  main="First order-absorption linear elimination", ## title for plot
  auto.key=list(adj=1,col=nlmixCOLS[c(2,1,3)],columns=3,space="bottom",rectangles=FALSE,points=FALSE) ## key for curves
)


#################################################################################
##                                                                             ##
## nlmixr2 analysis Part 2                                                      ##
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
run002F <-
  nlmixr2(One.comp.KA.ODE,          #the model definition
         PKdata,                   #the data set
         est = "focei",            #the estimation algorithm (FOCEi)
                                   #FOCEi options:
         foceiControl(print = 5))  #only print every 5th estimation step

## results are stored in the nlmixr object and can be viewed:
run002F

## and saved for future use or reference:
save(run002F, file = "run002F.Rdata")

ctr2 <- pmx_nlmixr(run002F, conts = c("WT","AGE"),cats=c("SEX","SPARSE"), vpc=FALSE,settings=pmx_settings(is.draft=FALSE))


ctr2 %>% pmx_plot_npde_pred
ctr2 %>% pmx_plot_npd_pred
ctr2 %>% pmx_plot_cwres_pred

ctr2 %>% pmx_report(name="ggPMX_report002F",
                   save_dir=".",
                   format="report",
                   extension="word")




#################################################################################
##                                                                             ##
##  Hands-on assignment: nlmixr2 model development                              ##
##  Examine the GOF plots and implement alternative absorption models          ##
##  -one or more transit compartment(s)                                        ##
##  - See slide for what a transit compartment is
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



