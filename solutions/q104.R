## Exercise 104
## Created: 2021-03-08 Enoch Chen 
## Edited:  2021-03-08 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q3.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr) # manipulate data

## Read data
## subset of stage 1
melanoma <- read_dta("melanoma.dta")%>% 
            filter(stage == 1)%>% 
            mutate(year8594 = ifelse(year8594 == 1, "Diagnosed 85-94", "Diagnosed 75-84"))
                   

#' stage (Clinical stage at diagnosis) is coded as follows
#' 0 [Unknown]
#' 1 [Localised]
#' 2 [Regional]
#' 3 [Distant]

#' status is coded as follows
#' 1 [Dead: cancer]
#' 2 [Dead: other]
#' 0 [Alive]
#' 4 [Lost to follow-up]


#' year8504 is coded as follows
#' 1 [Diagnosed 75-84]
#' 0 [Diagnosed 85-94]

##(a) Kaplan-Meier by calendar time 
mfityear8594 <- survfit(Surv(surv_mm, status==1) ~ year8594, data = melanoma)

plot(mfityear8594, col = 1:2,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
legend("bottomleft", legend=c(levels(as.factor(melanoma$year8594))), col=1:2, lty = 1)

##(b) Hazard function
plot(muhaz2(Surv(surv_mm,status==1)~year8594, data=melanoma))

##(c) Log-rank test 
survdiff(Surv(surv_mm, status==1) ~ year8594, data=melanoma)

## Wilcoxon test
survdiff(Surv(surv_mm, status==1) ~ year8594, data=melanoma, rho=1)

##(d) Estimate mortality rates

#' agegrp is coded as follows
#' 0 [0-44]
#' 1 [45-59]
#' 2 [60-74]
#' 3 [75+]

## Be aware of the unit of the rate: per 1000 person-months
survRate(Surv(surv_mm/1000, status==1)~agegrp, data=melanoma)

## Kaplan-Meier survival curves by agegrp
mfit_agegrp <- survfit(Surv(surv_mm, status==1) ~ agegrp, data = melanoma)

plot(mfit_agegrp, col = 1:4,
     xlab = "Months since diagnosis",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
legend("bottomright", legend=c(levels(as.factor(melanoma$agegrp))), col=1:4, lty = 1)

##(e) Re-scale time from months to years
mfit_agegrp_year <- survfit(Surv(surv_mm/12, status==1) ~ agegrp, data = melanoma)
plot(mfit_agegrp_year, col = 1:4,
     xlab = "Years since diagnosis",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
legend("bottomright", legend=c(levels(as.factor(melanoma$agegrp))), col=1:4, lty = 1)

## Be aware of the unit of the rate: per 1000 person-years
survRate(Surv(surv_mm/12/1000, status==1)~agegrp, data=melanoma)

##(f) Kaplan-Meier survival curves by sex

#' sex is coded as follows
#' 1 [Male]
#' 2 [Female]

mfit_sex <- survfit(Surv(surv_mm, status==1) ~ sex, data = melanoma)
plot(mfit_sex, col = 1:2,
     xlab = "Months since diagnosis",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates")
legend("bottomright", legend=c(levels(as.factor(melanoma$sex))), col=1:2, lty = 1)

## Hazard function
plot(muhaz2(Surv(surv_mm, status==1) ~ sex, data = melanoma))

## Estimate mortality rates by sex
## Be aware of the unit of the rate: per 1000 person-months
survRate(Surv(surv_mm/1000, status==1) ~ sex, data=melanoma)

## Log-rank test 
survdiff(Surv(surv_mm/1000, status==1) ~ sex, data=melanoma)
