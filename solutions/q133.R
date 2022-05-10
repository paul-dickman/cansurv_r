## Exercise 133
## Created: 2022-01-05 Enoch Chen 
## Edited:  2022-01-05 Enoch Chen
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr)   # manipulate data
library(rstpm2)
library(ggplot2) # ggplot

## Read data
melanoma <- read_dta("melanoma.dta")%>% 
            mutate(year8594 = ifelse(year8594 == 1, "Diagnosed 85-94", "Diagnosed 75-84"),
                   female = ifelse(sex == 2, +1, +0), # Create a female variable
                   event  = ifelse((status == 1) & surv_mm<= 60.5, 1, 0), # Redefine status
                   surv_mm = ifelse(surv_mm<= 60.5, surv_mm, 60.5),       # censored everyone after 60 months
                   t = surv_mm/12,                                     # surv_mm in years
                   agegrp1 = (agegrp== 0)+0, # used by time-dependent effect
                   agegrp2 = (agegrp== 1)+0, # used by time-dependent effect
                   agegrp3 = (agegrp== 2)+0, # used by time-dependent effect
                   agegrp4 = (agegrp== 3)+0) # used by time-dependent effect

melanoma <- melanoma %>% 
            mutate(year8594 = as_factor(year8594),
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

##(a) Proportional Hazards model
#' g(S)=log(-log(S))
ph <- stpm2(Surv(t, event) ~ female + year8594 + agegrp2 + agegrp3 + agegrp4, data = melanoma, df=4)
summary(ph)
eform(ph)

## (b) Proportional Odds Model
#' g(S)=-logit(S)
po <- stpm2(Surv(t, event) ~ female + year8594 + agegrp2 + agegrp3 + agegrp4, data = melanoma, df=4,
               link.type = "PO")
summary(po)
eform(po)

## (c) Compare survival and hazard function
## Prepare empty datasets for prediction
years <- levels(melanoma$year8594)
agegrps <- levels(melanoma$agegrp)

newdata1 <- data.frame(female=0, year8594=years[1], agegrp2=0, agegrp3=0, agegrp4=0)
newdata4 <- data.frame(female=0, year8594=years[1], agegrp2=0, agegrp3=0, agegrp4=1)

## Survival function
plot(ph, newdata=newdata1, 
     xlab="Time since diagnosis (years)",
     ylab="Survival probability",
     type="surv", ci=FALSE, ylim=c(0,1))
lines(po, newdata=newdata1, type="surv", col=2)
lines(ph, newdata=newdata4, type="surv", col=3)
lines(po, newdata=newdata4, type="surv", col=4)
legend("bottomright", legend = c("PH agegrp0", "PO agegrp0", "PH agegrp4", "PO agegrp4"), 
       col = 1:4, lty = 1)

## Hazard function
plot(ph, newdata=newdata1, 
     xlab="Time since diagnosis (years)",
     ylab="Hazard rates",
     type="haz", ci=FALSE, ylim=c(0,0.3))
lines(po, newdata=newdata1, type="haz", col=2)
lines(ph, newdata=newdata4, type="haz", col=3)
lines(po, newdata=newdata4, type="haz", col=4)
legend("topright", legend = c("PH agegrp0", "PO agegrp0", "PH agegrp4", "PO agegrp4"), 
                   col = 1:4, lty = 1)

##(d) Compare AIC and BIC 
AIC(ph, po)
BIC(ph, po)

##(e) Hazard ratio for female
plot(po, newdata=newdata1, 
     xlab="Time since diagnosis (years)",
     ylab="Hazard ratio",
     var="female",
     type="hr", ci=TRUE, ylim=c(0.45,0.65))

##(f) Compare hazard ratios for different covariate patterns	
plot(po, newdata=newdata1, 
     xlab="Time since diagnosis (years)",
     ylab="Hazard ratio",
     var="female",
     type="hr", ci=FALSE, ylim=c(0.5,0.7))
lines(po, newdata=newdata4, 
      var="female",
      type="hr", ci=FALSE, col=2)

##(g) Fit Aranda-Ordaz link function
ao <- stpm2(Surv(t, event) ~ female + year8594 + agegrp2 + agegrp3 + agegrp4, data = melanoma, df=4,
            link.type = "AO", theta.AO = 0)
summary(ao)
eform(ao)

AIC(ph, po, ao)
BIC(ph, po, ao)

##(h) Show estimate of theta with 95% CI
#' "spec" needs to be specified
lincom(po, "spec" ,eform = TRUE)

#' Code ends here


