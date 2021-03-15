## Exercise 111
## Created: 2021-03-15 Enoch Chen 
## Edited:  2021-03-15 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q6.html
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

#' status is coded as follows
#' 1 [Dead: cancer]
#' 2 [Dead: other]
#' 0 [Alive]
#' 4 [Lost to follow-up]

##(a)-i. Plot Kaplan-Meier curve
kaplanmeier <- survfit(Surv(surv_mm, event = status == 1) ~ year8594,
               data = melanoma)

plot(kaplanmeier,
     mark.time=F,
     xscale=12, # Turn months into years
     xlab="Years since diagnosis",
     ylab="S(t)",
     main = "Kaplan-Meier survival estimates",
     col=c("blue","red"))
# legend
legend("bottomleft",
       legend=c(levels(as.factor(melanoma$year8594))),
       col=c("blue","red"), lty = 1)

##(a)-ii. Hazard function
hazard <- muhaz2(Surv(surv_mm/12, status == 1) ~ year8594, data=melanoma)

plot(hazard,
     xlab="Years since diagnosis",
     ylab="Per person-years",
     main = "Smoothed hazard estimates",
     col=c("blue","red"),
     legend = FALSE)

legend("topright",
       legend=c(levels(as.factor(melanoma$year8594))),
       col=c("blue","red"), lty = 1)

##(a)-iii. Plot K-M curve and hazard function in one plot
par(mfrow=c(1,2)) ## Two graphs in the same window

# K-M
plot(kaplanmeier,
     mark.time=F,
     xscale=12, # Turn months into years
     xlab="Years since diagnosis",
     ylab="S(t)",
     main = "Kaplan-Meier survival estimates",
     col=c("blue","red"))
# legend
legend("bottomleft",
       legend=c(levels(as.factor(melanoma$year8594))),
       col=c("blue","red"), lty = 1)

# Hazard
plot(hazard,
     xlab="Years since diagnosis",
     ylab="Per person-years",
     main = "Smoothed hazard estimates",
     col=c("blue","red"),
     legend = FALSE)

legend("topright",
       legend=c(levels(as.factor(melanoma$year8594))),
       col=c("blue","red"), lty = 1)

##(b) Hazard rate and Poisson regression
# Per person-years
survRate(Surv(surv_mm/12/1000,  status == 1) ~ year8594, data=melanoma)

poisson_b <- glm( status == 1 ~ year8594 + offset( log( surv_mm/12/1000 ) ), family=poisson, data=melanoma)
summary(poisson_b)
eform(poisson_b)

##(c)-i. Restrict follow-up to 10 years
melanoma_10y <- mutate(melanoma,
                       status= ifelse(surv_mm <= 120 & status == 1, 1, 0),
                       surv_mm = ifelse(surv_mm<=120, surv_mm, 120) # censored everyone after 120 months
)

##(c)-ii.  Hazard rate and Poisson regression
survRate(Surv(surv_mm/12/1000, status == 1) ~ year8594, data=melanoma_10y)

# Calculate by hand the ratio (85-94 to 75-84)
27.778/31.453 # 0.883159

##(c)-iii. Poisson regression
poisson_c <- glm( status == 1 ~ year8594 + offset( log( surv_mm/12/1000  ) ), family=poisson, data=melanoma_10y)
summary(poisson_c)
eform(poisson_c)

