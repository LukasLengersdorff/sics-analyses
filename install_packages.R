
#############################################################
### Install necessary packages ##############################
#############################################################

## --- Prototype SICS package ---
## 
## This is used for analyses, and installed from Github

# devtools is needed to install from Github
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# SICS
devtools::install_github("LukasLengersdorff/sics-prototype")


## --- Other packages ---
## 

packages = c("tidyverse", # for convenience
             "osfr",      # to load example data from OSF
             "mvtnorm",   # used in sampling demo for MV normal density
             "ggpubr","gridExtra", # Used to arrange several ggplots
             "knitr")     # For LateX tables


new_packages = packages[!packages %in% installed.packages()[, "Package"]]

if (length(new_packages)) {
  install.packages(new_packages)
}