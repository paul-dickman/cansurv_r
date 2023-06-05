## Exercise 103
## Created: 2021-03-04 Enoch Chen 
## Edited:  2023-06-05 Enoch Chen: add "list the first few observations"
##          2021-03-04 Enoch Chen 
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q2.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)

## Read data
melanoma <- read_dta("melanoma.dta")

#' status is coded as follows
#' 1 [Dead: cancer]
#' 2 [Dead: other]
#' 0 [Alive]
#' 4 [Lost to follow-up]

#' stage (Clinical stage at diagnosis) is coded as follows
#' 0 [Unknown]
#' 1 [Localised]
#' 2 [Regional]
#' 3 [Distant]

## List the first few observations
head(melanoma)

## Tabulate stage distribution
Freq <- xtabs(~stage, data=melanoma)
cbind(Freq, Percent=prop.table(Freq) %>% `*`(100) %>% round(2))

##(a) Survival and hazard function by stage
## Survival
par(mfrow=c(1, 2)) # Will plot 2 graphs
mfit <- survfit(Surv(surv_mm, status == 1) ~ stage, data = melanoma)

plot(mfit, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival")

## Hazard
hazards <- muhaz2(Surv(surv_mm, status == 1)~stage, melanoma)
plot(hazards,
     col=1:4, lty=1, xlim=c(0,250), ylim=c(0,0.08),
     legend.args=list(bty="n"))


##(b) Mortality rates by stage 
## events/person-months (be aware of the unit!)
survRate(Surv(surv_mm, status == 1) ~ stage, data = melanoma)

## events/person-years 
survRate(Surv(surv_mm/12, status == 1) ~ stage, data = melanoma)


##(c) events/1000 person-years
survRate(Surv(surv_mm/12/1000, status == 1) ~ stage, data = melanoma)


##(d)  Survival and hazard function by sex

## sex was coded as follows
#' 1 [Male]
#' 2 [Female]

## Survival
par(mfrow=c(1, 2)) # Will plot 2 graphs
mfit <- survfit(Surv(surv_mm, status == 1) ~ sex, data = melanoma)

plot(mfit, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival")

## Hazard
hazards <- muhaz2(Surv(surv_mm, status == 1)~ sex, melanoma)
plot(hazards,
     col=1:2, lty=1)

## Mortality rates (per 1000 person-years) by sex 
survRate(Surv(surv_mm/12/1000, status == 1) ~ sex, data=melanoma)

##(e) Tabulate status by age group

## agegrp is coded as follows
#' 0 [0-44]
#' 1 [45-59]
#' 2 [60-74]
#' 3 [75+]

xtabs(~status+agegrp, melanoma)

##(f) All-cause survival by stage
par(mfrow=c(1, 1))
afit <- survfit(Surv(surv_mm, status == 1|status == 2) ~ stage, data = melanoma)
plot(afit, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates\nAll-cause")
legend("topright", legend=c(levels(as.factor(melanoma$stage))), col=1:4, lty = 1)

##(g) Compare cause-specific and all cause survival
## agegrp==3 is age 75+
par(mfrow=c(1, 2))
mfit75 <- survfit(Surv(surv_mm, status == 1) ~ stage, data = subset(melanoma,agegrp==3))

plot(mfit75, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates\nCancer | Age 75+")
legend("topright", legend=c(levels(as.factor(melanoma$stage))), col=1:4, lty = 1)

afit75 <- survfit(Surv(surv_mm, status == 1|status == 2) ~ stage, data = subset(melanoma,agegrp==3))

plot(afit75, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier survival estimates\nAll-cause | Age 75+")
legend("topright", legend=c(levels(as.factor(melanoma$stage))), col=1:4, lty = 1)

##(h) Estimate both cancer-specific and all-cause mortality for each age group 
par(mfrow=c(1, 2))
mfitage <- survfit(Surv(surv_mm, status == 1) ~ agegrp, data = melanoma)
plot(mfitage, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier estimates of\ncancer survival by age group")
legend("bottomright", legend=c(levels(as.factor(melanoma$agegrp))), col=1:4, lty = 1)

afitage <- survfit(Surv(surv_mm, status == 1|status == 2) ~ agegrp, data = melanoma)
plot(afitage, col=1:4,
     xlab = "Follow-up Time",
     ylab = "Survival",
     main = "Kaplan-Meier estimates of\nall-cause survival by age group")
legend("topright", legend=c(levels(as.factor(melanoma$agegrp))), col=1:4, lty = 1)