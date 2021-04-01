## Exercise 123
## Created: 2021-04-01 Enoch Chen 
## Edited:  2021-04-01 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q12.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr)   # manipulate data

## Read data
## All stages
melanoma <- read_dta("melanoma.dta")%>% 
            mutate(year8594 = ifelse(year8594 == 1, "Diagnosed 85-94", "Diagnosed 75-84"),
                   agegrp = as.factor(agegrp),
                   sex = as.factor(sex),
                   stage = as.factor(stage),
                   subsite = as.factor(subsite))

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
#' 
#' subsite is coded as follows
#' 1  [Head and Neck]
#' 2  [Trunk]
#' 3  [Limbs]
#' 4  [Multiple and NOS]
#' 
#' stage is coded as follows
#' 0  [Unknown]
#' 1  [Localised]
#' 2  [Regional] 
#' 3  [Distant]

# Create new variable for cause-specific death
melanoma <- mutate(melanoma,
                       status_cancer = ifelse(status == 1 , 1, 0),
                       )

##(a) Cox regression
# Note:  R uses the Efron method for approximating the likelihood in the presence of ties.
# Stata uses the Breslow method.
cox_a <- coxph(Surv(surv_mm, status_cancer) ~ sex, data=melanoma)
summary(cox_a)


##(b) Cox regression controlling for possible confounders
cox_b <- coxph(Surv(surv_mm, status_cancer) ~ sex + agegrp + stage + subsite + year8594, data=melanoma)
summary(cox_b)


##(c) Cox regression with interactions
cox_c <- coxph(Surv(surv_mm, status_cancer) ~ sex * agegrp, data=melanoma)
summary(cox_c)

# Test of whether the HR's for sex differ across age groups
anova(cox_c)

##Cox regression with main effects and interactions
cox_c2 <- coxph(Surv(surv_mm, status_cancer) ~ year8594 + subsite + stage + sex * agegrp, data=melanoma)
summary(cox_c2)

# Test of whether the HR's for sex differ across age groups
anova(cox_c2)


##(d) Cox regression, "best model"
cox_d <- coxph(Surv(surv_mm, status_cancer) ~ sex + year8594 + agegrp + subsite + stage, data=melanoma)
summary(cox_d)
cox.zph(cox_d, transform="identity") # Stata default

# Stratification
cox_strat1 <- coxph(Surv(surv_mm, status_cancer) ~ sex + year8594 + agegrp + subsite + strata(stage), data=melanoma)
summary(cox_strat)
cox.zph(cox_strat1, transform="identity") # Stata default

cox_strat2 <- coxph(Surv(surv_mm, status_cancer) ~ sex + year8594 + agegrp  + strata(stage), data=melanoma)
summary(cox_strat2)
cox.zph(cox_strat2, transform="identity")

# Additional interaction between age at diagnosis and stage.
cox_strat3 <- coxph(Surv(surv_mm, status_cancer) ~ sex * agegrp + year8594 + agegrp + subsite + strata(stage), data=melanoma)
summary(cox_strat3)
anova(cox_strat3)
