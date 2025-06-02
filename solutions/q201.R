## Exercise 201
## Created: 2023-06-05 Enoch Chen, Joshua Entrop
## Edited:  2023-06-06 Enoch Chen: Debug status
## Updated: 2024-11-18 Paul Dickman (add script to create ratetable)
###############################################################################
# Load packages
library(haven)
library(dplyr)
library(tidyr)
library(lubridate)
library(relsurv)

# Load melanoma data for limited stages
melanoma <- read_dta("../data/melanoma.dta") %>% 
            filter(stage == 1) # Localised

# Prepare melanoma dataset
melanoma <- melanoma %>% 
            mutate(surv_dd = surv_mm * (365.241/12), # Time in days for relsurv pkg
                   status  = if_else(status %in% c(1, 2), 1, 0),
                   year    = year(dx))

# Create an R ratetable object
source("make_ratetable.R")

# Estimate relative survival
# (a) Annual intervals
rs_fit_e2 <- rs.surv(Surv(surv_dd, status == 1 )~ year8594,
                     data      = melanoma,
                     ratetable = ratetable,
                     method    = "ederer2",
                     precision = 12,
                     rmap      = list(age = age * 365.241))

summary(rs_fit_e2, times = c(5, 10) * 365.241)

# (b) 6 month intervals
rs_fit_e2 <- rs.surv(Surv(surv_dd, status == 1 )~ 1,
                     data      = melanoma,
                     ratetable = ratetable,
                     method    = "ederer2",
                     precision = 365.241/2,
                     rmap      = list(age = age * 365.241))

summary(rs_fit_e2, times = c(5, 10) * 365.241)

# (c) 3 month intervals
rs_fit_e2 <- rs.surv(Surv(surv_dd, status == 1 )~ 1,
                     data      = melanoma,
                     ratetable = ratetable,
                     method    = "ederer2",
                     precision = 365.24/4,
                     rmap      = list(age = age * 365.24))
summary(rs_fit_e2, times = c(5, 10) * 365.241)

# (d) annual intervals up to 20 years
rs_fit_e2 <- rs.surv(Surv(surv_dd, status == 1 )~ 1,
                     data      = melanoma,
                     ratetable = ratetable,
                     method    = "ederer2",
                     precision = 365.24,
                     rmap      = list(age = age * 365.24))
summary(rs_fit_e2, times = c(5, 10, 20) * 365.241)

# (e) Plot by period
plot(rs_fit_e2)

rs_fit_e2 <- rs.surv(Surv(surv_dd, status == 1 )~ year8594,
                     data      = melanoma,
                     ratetable = ratetable,
                     method    = "ederer2",
                     precision = 365.24,
                     rmap      = list(age = age * 365.24))

# (f) To be completed
# (g) (h)
rs_fit_e1 <- rs.surv(Surv(surv_dd, status == 1 )~ 1,
                     data      = melanoma,
                     ratetable = ratetable,
                     method    = "ederer1",
                     precision = 365.24,
                     rmap      = list(age = age * 365.24))
summary(rs_fit_e1, times = c(5, 10) * 365.241)

rs_fit_ha <- rs.surv(Surv(surv_dd, status == 1 )~ 1,
                   data      = melanoma,
                   ratetable = ratetable,
                   method    = "hakulinen",
                   precision = 365.24,
                   rmap      = list(age = age * 365.24))
summary(rs_fit_ha, times = c(5, 10) * 365.241)

rs_fit_pp <- rs.surv(Surv(surv_dd, status == 1 )~ 1,
                     data      = melanoma,
                     ratetable = ratetable,
                     method    = "pohar-perme",
                     precision = 365.24,
                     rmap      = list(age = age * 365.24))
summary(rs_fit_pp, times = c(5, 10) * 365.241)

# (i) To be completed
