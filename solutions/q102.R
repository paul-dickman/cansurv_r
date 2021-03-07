## Exercise 102
## Created: 2021-03-03 Enoch Chen 
## Edited:  2021-03-03 Enoch Chen 
###############################################################################
## Load the packages
library(biostat3)
library(dplyr)
library(haven)

## Read data
melanoma <- read_dta("melanoma.dta")

#' Vital status (status) is coded as follows
#' 1 [Dead: cancer]
#' 2 [Dead: other]
#' 0 [Alive]
#' 4 [Lost to follow-up]

## Generate a new failure variable
melanoma <- transform(melanoma,
                      csr_fail = ifelse( status == 1, 1, 0))

## Life table
print(lifetab2(Surv(floor(surv_yy), status == 1)~1, melanoma, breaks=0:20), digits=2)
print(lifetab2(Surv(floor(surv_mm), status == 1)~1, melanoma, breaks=c(12*(0:20))), digits=2)

## Kaplan-Meier (survival time in years)
mfit1 <- survfit(Surv(surv_yy, status == 1) ~ 1, data = melanoma)
summary(mfit1)          

## Kaplan-Meier (survival time in months)
mfit2 <- survfit(Surv(surv_mm, status == 1) ~ 1, data = melanoma)
summary(mfit2)          