##------------------ Introduction --------------------------
#
# The purpose of this document is to to generate the diagnostic plots using ggPMX for the warfarin 
# PKPD model using nlmixr. 
# The material is designed to run on R 4.1.0 with nlmixr2 v2.0.6 and ggPMX v1.2.8.
#
#------------------------------------------------------------

rm(list=ls(all=TRUE)) # Clean up environment

# set-up the model output directory - the folder containing modeling outputs
work_dir = "~/ggPMX/trunk/PUBLICATION/PAGE2022"
knitr::opts_chunk$set(echo=FALSE,root.dir=work_dir)

# set-up working directory
setwd(work_dir)

# To be executed only once for a given R session
if(FALSE){ # press enter when asked a question
  withr::with_libpaths("~/ggPMX/trunk/PUBLICATION/PAGE2022/libs/",{
    remotes::install_github("nlmixr2/lotri")
    remotes::install_github("nlmixr2/rxode2")
    remotes::install_github("nlmixr2/nlmixr2data")
    remotes::install_github("nlmixr2/nlmixr2est")
    remotes::install_github("nlmixr2/nlmixr2extra")
    remotes::install_github("nlmixr2/nlmixr2plot")
    remotes::install_github("nlmixr2/nlmixr2")
    remotes::install_github("ggPMXdevelopment/ggPMX")
  }, action = "prefix")
}

withr::with_libpaths("~/ggPMX/trunk/PUBLICATION/PAGE2022/libs/",{
  library(rlang)
  library(rxode2)
  library(nlmixr2est)
  library(nlmixr2extra)
  library(nlmixr2plot)
  library(nlmixr2)
  library(ggPMX)
}, action = "prefix")

library(dplyr)
library(ggplot2)

# Check library versions
sessionInfo()

#------------------------------------------------------------
## Model Fitting
# 
# **Exercise 1:** Create a PKPD model using nlmixr2 (use the demo example below):
# 
# * Create the warfarin PKPD model using nlmixr2 (use demo code below)
# * Fit the model to warfarin data (warfarin dataset is included by default in the nlmixr2 installation; it contains weight (wt), age and sex as covariates)
# * Save the model fit as fit.rds in the working directory
# 
# Alternative: you can also use Monolix or NONMEM with any model if you wish. To best explore ggPMX capacities, it would be preferable if there are continuous and categorical covariates in your dataset.

# SOLUTION

# Warfarin example
pk.turnover.emax3 <- function() {
  ini({
    tktr <- log(1)
    tka <- log(1)
    tcl <- log(0.1)
    tv <- log(10)
    ##
    eta.ktr ~ 1
    eta.ka ~ 1
    eta.cl ~ 2
    eta.v ~ 1
    prop.err <- 0.1
    pkadd.err <- 0.1
    ##
    temax <- logit(0.8)
    tec50 <- log(0.5)
    tkout <- log(0.05)
    te0 <- log(100)
    ##
    eta.emax ~ .5
    eta.ec50  ~ .5
    eta.kout ~ .5
    eta.e0 ~ .5
    ##
    pdadd.err <- 10
  })
  model({
    ktr <- exp(tktr + eta.ktr)
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    emax = expit(temax+eta.emax)
    ec50 =  exp(tec50 + eta.ec50)
    kout = exp(tkout + eta.kout)
    e0 = exp(te0 + eta.e0)
    ##
    DCP = center/v
    PD=1-emax*DCP/(ec50+DCP)
    ##
    effect(0) = e0
    kin = e0*kout
    ##
    d/dt(depot) = -ktr * depot
    d/dt(gut) =  ktr * depot -ka * gut
    d/dt(center) =  ka * gut - cl / v * center
    d/dt(effect) = kin*PD -kout*effect
    ##
    cp = center / v
    cp ~ prop(prop.err) + add(pkadd.err)
    effect ~ add(pdadd.err) | pca
  })
}
ui3 = nlmixr(pk.turnover.emax3)
ui3

# Dataset used to fit the model (imbedded in nlmixr2):
summary(warfarin)

# Fit the model to the data
fit <- nlmixr(pk.turnover.emax3, warfarin, "saem", control=list(print=0),
              table=list(cwres=TRUE, npde=TRUE))
print(fit)

saveRDS(fit,"fit.rds")
#------------------------------------------------------------------

## ggPMX Controller
# 
# **Exercise 2:** ggPMX controller creation
# 
# * Create the controller **for PK only (DVID=="cp")**
# * Plot NPDE vs TIME
# * Plot NPDE vs TIME for time > 50
# * Plot NPDE vs PRED and change the labels of x and y axes
# * Plot DV vs PRED on log-log scale
# * Plot DV vs IPRED and stratify it by categorical covariates
# * Plot EBE vs continuous covariates and EBE vs categorical covariates
# * Stratify NPDE vs TIME by categorical covariates
# 
# Hints:
# 
# * To get help use `?pmx_plot_npde_pred`, `?pmx_plot_npde_time`, etc.

# SOLUTION

# Load the fit
fit <- readRDS("fit.rds")
capture.output(print(fit)) %>%
  sapply(crayon::strip_style) %>%
  unname %>%
  print

head(warfarin) # check columns of dataset
unique(warfarin$dvid) # check values of dvid

# Define the endpoint
my_endpoint = pmx_endpoint(code = "cp")

# Create controller for PK (default):
ctr = pmx_nlmixr(fit, endpoint = my_endpoint, conts = c("wt", "age"), cats = "sex")

# Check what is in ctr
ctr

# Plot NPDE vs TIME:
ctr %>% pmx_plot_npde_time()

# Plot NPDE vs TIME for time > 50
ctr %>% pmx_plot_npde_time(filter = time > 50)

# Plot NPDE vs PRED and change the labels of x and y axes
ctr %>% pmx_plot_npde_pred(labels = list(x = "My new x-label", y = "My new y-label"))

# Plot DV vs PRED on log-log scale:
ctr %>% pmx_plot_dv_pred
ctr %>% pmx_plot_dv_pred(scale_x_log10 = TRUE, scale_y_log10 = TRUE)

# Plot DV vs IPRED on log-log scale, stratified by categorical covariates
ctr %>% pmx_plot_dv_ipred(scale_x_log10 = TRUE, scale_y_log10 = TRUE, strat.facet = "sex")

# Stratify NPD vs TIME colored by age
ctr %>% pmx_plot_npd_time(strat.color = "age")

# Plot EBE vs continuous covariates and EBE vs categorical covariates
ctr %>% pmx_plot_eta_conts()
ctr %>% pmx_plot_eta_cats()
#---------------------------------------------------------------------------

## ggPMX Report
# 
# **Exercise 3:** ggPMX report generation - **for PK only**
# 
# * Create a folder named GOF_PK where will be stored the goodness-of-fit plots for the PK
# * Generate Word and pdf ggPMX reports (named *ggPMX_report_pk.docx* and *ggPMX_report_pk.pdf*) in the GOF_PK folder; use the option generating the diagnostic plots folder that contains each figure as a separate file.
# * Generate an HTML report without footnote (the footnote indicates the source location of each figure) and without the folder with figures as separate files.
# 
# Hint: To get help on the various options, use `?pmx_report`.

# SOLUTION

# File Name for ggPMX report
report_name = "ggPMX_report_pk"

my_save_dir = "GOF_PK"
system(paste0("mkdir ", my_save_dir))

# Generate the pdf and docx ggPMX report, as well as plots folder
ctr %>% pmx_report(name      = report_name, 
                   save_dir  = file.path(work_dir, my_save_dir), 
                   format    = "both", 
                   footnote  = TRUE, 
                   extension = c("word"))

# Generate the html ggPMX report WITHOUT GoF plots folder and NO footnote
ctr %>% pmx_report(name      = report_name, 
                   save_dir  = file.path(work_dir, my_save_dir), 
                   format    = "report", 
                   footnote  = FALSE, 
                   extension ="html")
#-----------------------------------------------------------------
# 
# **Exercise 4:** ggPMX Report Customization
# 
# * Edit your own template based on the file *ggPMX_report_pk.Rmd*.
# * Set format to landscape in the new template.
# * Change the color of the loess line to blue for the NPDE vs PRED plot and change the axis labels.
# * Remove the loess line in the IWRES vs TIME plot.
# * Filter TIME>30 in NPDE vs TIME plot.
# * Set log-x and log-y axes for the DV vs PRED plot.
# * Include all pages of individual plots (grid 4x4).
# * Save the newly created Rmarkown file in another folder (name it *ggPMX_report_pk_custom.Rmd*) and execute it using `pmx_report()`.
#  
# Hints:
# 
# * To get help use `?pmx_plot_npde_pred`, `?pmx_plot_npde_time`, `?pmx_plot_iwres_time`, `?pmx_plot_dv_pred` and `?pmx_plot_individual`.
# * Landscape format is set-up with `classoption: landscape` in the RMarkdown header

# SOLUTION

# File Name for ggPMX report
report_name_custom = "ggPMX_report_pk_custom"

# Add classoption: landscape in Rmd Header
# ---
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output: 
#   pdf_document:
#     toc: true
#   word_document:
#     toc: true
#   html_document:
#     toc: true
# classoption: landscape
# ---

# Edit the ggPMX report on RStudio editor and overwrite the respective chunks with the following:
ctr %>% pmx_plot_npde_pred(smooth = list(color='blue'), labels = list(x = "My PRED", y = "My NPDE"))
ctr %>% pmx_plot_npde_time(filter = TIME > 30)
ctr %>% pmx_plot_iwres_time(is.smooth = FALSE)
ctr %>% pmx_plot_dv_pred(scale_x_log10 = TRUE, scale_y_log10 = TRUE)
ctr %>% pmx_plot_individual(facets = list(nrow = 4, ncol = 4), which_pages = "all")

# Define the name of the custom report file
report_name_custom = "ggPMX_report_pk_custom"
ctr %>% pmx_report(name      = report_name_custom, 
                   template  = "ggPMX_report_pk_custom.Rmd", # has to be saved in different directory
                   save_dir  = file.path(work_dir, my_save_dir), 
                   format    = "report",
                   footnote  = TRUE, 
                   extension = "all")
#---------------------------------------------------------------------
# 
# **Exercise 6:** ggPMX Report for PD
# 
# * Create a folder named GOF_PD where will be stored the goodness-of-fit plots for the PD
# * Generate a Word ggPMX report (named *ggPMX_report_pd.docx*) and the diagnostic plots folder for the PD endpoint (dvid=="pca")
# * Remove DRAFT label on all plots
# 
# Hints:
# 
# * First need to create the controller for PD endpoint (dvid=="pca")
# * Generate the report using the PD controller with the command `pmx_report()`; output files should be stored in folder GOF_PD
# * Explore `?pmx_settings`

# SOLUTION

my_save_dir_pd = "GOF_PD"
system(paste0("mkdir ", my_save_dir_pd))

# File Name for ggPMX PD report
report_name_pd = "ggPMX_report_pd"

# Create the controller for PD endpoint (DVID==2)
# Draft setting (TRUE or FALSE) performed when creating the controller
ctrPD = pmx_nlmixr(fit, 
                   endpoint = pmx_endpoint(code = "pca"),
                   conts = c("wt", "age"), 
                   cats = "sex",
                   settings = pmx_settings(is.draft=FALSE))

# Generate the pdf ggPMX report and diagnostic plots folder
ctrPD %>% pmx_report(name       = report_name_pd, 
                     save_dir   = file.path(work_dir, my_save_dir_pd), 
                     format     = "both",
                     extension  = "word")

#----------------------- END ------------------------------------------------
