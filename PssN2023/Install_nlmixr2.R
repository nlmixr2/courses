# For Windows, you must install Rtools first
# Then install the following packages with their dependencies:
install.packages(c("nlmixr2","ggPMX","xpose.nlmixr2","patchwork"),dependencies=TRUE)

##Check to see if nlmixr2 is installed correctly:
library(nlmixr2)
nlmixr2CheckInstall()

one.compartment.IV.model <- function(){
  ini({
    lCl <- 1.6      #log Cl (L/hr)   
    lVc <- 4.5      #log V (L)   
    prop.err <- 0.3 
    eta.Vc ~ 0.1   
    eta.Cl ~ 0.1   
  })
  model({
    Vc <- exp(lVc + eta.Vc)
    Cl <- exp(lCl + eta.Cl)
    d / dt(centr) = -(Cl / Vc) * centr;
    cp = centr / Vc;
    cp ~ prop(prop.err)
  })
}

fit <- nlmixr2(one.compartment.IV.model, Bolus_1CPT,  est="focei")

fit$OBJF
