## Exercise 200
## Created: 2023-06-05 Enoch Chen 
## Edited:  2023-06-05 Enoch Chen
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(relsurv)
library(tidyverse)

## Life table
colon_sample <- biostat3::colon_sample


# Set up the data for relative survival analysis

# Fit the relative survival model
rs_fit <- relsurv(surv_obj, popmort)
