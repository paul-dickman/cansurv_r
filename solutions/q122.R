## Exercise 122
## Created: 2021-04-01 Enoch Chen 
## Edited:  2021-04-01 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q11.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr)   # manipulate data

## Read data
## subset of stage 1
melanoma <- read_dta("melanoma.dta")%>% 
            filter(stage == 1)%>% 
            mutate(year8594 = ifelse(year8594 == 1, "Diagnosed 85-94", "Diagnosed 75-84"),
                   agegrp = as.factor(agegrp),
                   sex = as.factor(sex))

#' status is coded as follows
#' 1 [Dead: cancer]
#' 2 [Dead: other]
#' 0 [Alive]
#' 4 [Lost to follow-up]

#' agegrp is coded as follows
#' 1 [0-44]
#' 2 [45-59]
#' 3 [60-74]
#' 4 [75+]

#' sex is coded as follows
#' 1 [Male]
#' 2 [Female]

# Restrict follow-up to 10 years and create new variable for all-cause death and cause-specific death
melanoma_10y <- mutate(melanoma,
                       status_all= ifelse(surv_mm <= 120 & (status == 1 | status == 2), 1, 0),
                       status_cancer = ifelse(surv_mm <= 120 & status == 1 , 1, 0),
                       surv_mm = ifelse(surv_mm<=120, surv_mm, 120), # censored everyone after 120 months
                       )

##(a) Cox regression with all-cause death
cox_a <- coxph(Surv(surv_mm, status_all) ~ sex + year8594 + agegrp,
                   data = melanoma_10y)
summary(cox_a)

##(b) Cox regression with cause-specific death (cancer)
cox_b <- coxph(Surv(surv_mm, status_cancer) ~ sex + year8594 + agegrp,
               data = melanoma_10y)
summary(cox_b)