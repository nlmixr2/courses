# Setting up nlmixr2

# Install nlmixr2-verse 
# you get all you need to set up nlmixr2

install.packages("nlmixr2",dependencies = TRUE)

# Installation of supportive scripts

install.packages(c("xpose.nlmixr2", # Additional goodness of fit plots
                   # baesd on xpose
                   "nlmixr2targets", # Simplify work with the
                   # `targets` package
                   "babelmixr2", # Convert/run from nlmixr2-based
                   # models to NONMEM, Monolix, and
                   # initialize models with PKNCA
                   "nonmem2rx", # Convert from NONMEM to
                   # rxode2/nlmixr2-based models
                   "nlmixr2lib", # a model library and model
                   # modification functions that
                   # complement model piping
                   "nlmixr2rpt" # Automated Microsoft Word and
                   # PowerPoint reporting for nlmixr2
),dependencies = TRUE)



# check installation

nlmixr2::nlmixr2CheckInstall()

library(nlmixr2)
nlmixr2update()


# install additional handy packages

install.packages("tidyverse",dependencies = TRUE)
# install packages for survival analysis
install.packages("survminer",dependencies = TRUE)