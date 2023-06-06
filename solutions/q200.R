## Exercise 200
## Created: 2023-06-05 Enoch Chen 
## Edited:  2023-06-06 Enoch Chen
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(relsurv)
library(tidyverse)

## Life table
colon_sample <- biostat3::colon_sample

## Make ratetable
source("../solutions/make_ratetable.R")

# Set up the data for relative survival analysis
colon_sample <- colon_sample %>% 
                mutate(surv_dd = surv_mm * 30, # Time in days for relsurv pkg
                       status  = if_else(status %in% c("Dead: cancer", "Dead: other"), 1, 0))

rs_fit_e1 <- rs.surv(Surv(surv_dd, status == 1)~ 1,
                     data      = colon_sample,
                     ratetable = ratetable,
                     method    = "ederer1",
                     rmap      = list(age = age * 365.241, year=yydx))

rs_fit_e2 <- rs.surv(Surv(surv_dd, status == 1)~ 1,
                     data      = colon_sample,
                     ratetable = ratetable,
                     method    = "ederer2",
                     rmap      = list(age = age * 365.241, year=yydx))
## Caution! The results looked weird to me
## I suspected the estimates were wrong due to either ratetable or time in days?
summary(rs_fit_e1, times = c(0, 1, 2, 3, 4, 5) * 365.241)
summary(rs_fit_e2, times = c(0, 1, 2, 3, 4, 5) * 365.241)

