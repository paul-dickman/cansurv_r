## Exercise 132
## Created: 2021-08-05 Enoch Chen 
## Edited:  2021-08-05 Enoch Chen
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr)   # manipulate data
library(rstpm2)
library(ggplot2) # ggplot

## Read data
melanoma <- read_dta("melanoma.dta")%>% 
            filter(stage == 1) %>%
            mutate(year8594 = ifelse(year8594 == 1, "Diagnosed 85-94", "Diagnosed 75-84"),
                   female = ifelse(sex == 2, 1, 0), # Create a female variable
                   event  = ifelse((status == 1) & surv_mm<= 60.5, 1, 0), # Redefine status
                   surv_mm = ifelse(surv_mm<= 60.5, surv_mm, 60.5),       # censored everyone after 60 months
                   t = surv_mm/12,                                     # surv_mm in years
                   agegrp1 = (agegrp== 0)+0, # used by time-dependent effect
                   agegrp2 = (agegrp== 1)+0, # used by time-dependent effect
                   agegrp3 = (agegrp== 2)+0, # used by time-dependent effect
                   agegrp4 = (agegrp== 3)+0) # used by time-dependent effect

melanoma <- melanoma %>% 
            mutate(year8594 = as_factor(year8594),
                   female = as_factor(female),
                   agegrp = as_factor(agegrp)) 
str(melanoma)

#' status is coded as follows
#' 1 [Dead: cancer]
#' 2 [Dead: other]
#' 0 [Alive]
#' 4 [Lost to follow-up]

#' agegrp is coded as follows
#' 0 [0-44]  -> agegrp1
#' 1 [45-59] -> agegrp2
#' 2 [60-74] -> agegrp3
#' 3 [75+]   -> agegrp4

#' female is coded as follows
#' 0 [Male]
#' 1 [Female]

##(a) fit a Cox model 
# Assess the PH assumption for age group
cox_a <- coxph(Surv(t, event) ~ female + year8594 + agegrp2 + agegrp3 + agegrp4, data = melanoma)
summary(cox_a)

##Plot of the scaled Schoenfeld residuals for calendar period 1985–94
cox_a.phtest <- cox.zph(cox_a, transform="identity") 
print(cox_a.phtest)
plot(cox_a.phtest[3], resid=TRUE, se=TRUE, main="Test of PH assumption: Schoenfeld residuals")
plot(cox_a.phtest[4], resid=TRUE, se=TRUE, main="Test of PH assumption: Schoenfeld residuals")
plot(cox_a.phtest[5], resid=TRUE, se=TRUE, main="Test of PH assumption: Schoenfeld residuals")

##(b) fit flexible parametric model
fpm_b <- stpm2(Surv(t, event) ~ female + year8594 + agegrp2 + agegrp3 + agegrp4, data = melanoma, df=4)
summary(fpm_b)
eform(fpm_b)

## Hazard function
females <- levels(melanoma$female)
years <- levels(melanoma$year8594)
agegrps <- levels(melanoma$agegrp)

newdata1 <- data.frame(female=females[1], year8594=years[1], agegrp2=0, agegrp3=0, agegrp4=0)
newdata2 <- data.frame(female=females[1], year8594=years[1], agegrp2=1, agegrp3=0, agegrp4=0)
newdata3 <- data.frame(female=females[1], year8594=years[1], agegrp2=0, agegrp3=1, agegrp4=0)
newdata4 <- data.frame(female=females[1], year8594=years[1], agegrp2=0, agegrp3=0, agegrp4=1)

plot(fpm_b, newdata=newdata1, 
     xlab="Time since diagnosis (years)",
     ylab="Cause specific mortality rate (per pys)",
     type="haz", ci=FALSE, ylim=c(0,0.2))

lines(fpm_b, newdata=newdata2, type="haz", lty=2)
lines(fpm_b, newdata=newdata3, type="haz", lty=3)
lines(fpm_b, newdata=newdata4, type="haz", lty=4)
legend("topright", legend=paste0(agegrps), lty=1:4)
