## Exercise 111
## Created: 2021-03-15 Enoch Chen 
## Edited:  2021-03-15 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q7.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr) # manipulate data

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

##(d) Split follow-up into 1-year intervals
# Create a new time variable in years using surv_mm
melanoma_10y <- mutate( melanoma_10y, surv_yy1 = surv_mm/12)

# Split the survival time by year
melanoma_10y_spl <- survSplit(melanoma_10y, cut=0:9, end="surv_yy1", start="start",
                                event = "status" ) # event must be a charater variable

## Calculate follow-up time and recode start time as a factor
melanoma_10y_spl <- mutate(melanoma_10y_spl,
                           pt = surv_yy1 - start,
                           fu = as.factor(start) )

##(e) Tabulate (and produce a graph of) the rates by follow-up time 
# Calculate the hazard rate by observation year
yearly_rates <- survRate(Surv(pt/1000, status == 1) ~ fu, data=melanoma_10y_spl)

# Plot
with(yearly_rates, matplot(as.numeric(as.character(fu))+0.5,
                           cbind(rate, lower, upper),
                           ylim=c(0,max(upper)),
                           lty=c("solid","dashed","dashed"),
                           col=c("black","gray","gray"),
                           type="l",
                           main="Cancer deaths by years since diagnosis",
                           ylab="Hazard rate per 1000 person-years",
                           xlab="Years since diagnosis") )

##(f) Compare the plot of the estimated rates to a plot of the hazard rate
par(mfrow=c(1,2))
with(yearly_rates, matplot(as.numeric(as.character(fu))+0.5,
                           cbind(rate, lower,
                                 upper),
                           ylim=c(0,max(upper)),
                           lty=c("solid","dashed","dashed"),
                           col=c("black","gray","gray"),
                           type="l",
                           main="Cancer deaths by time since diagnosis",
                           ylab="Hazard rate per 1000 person-years",
                           xlab="Years since diagnosis") )

hazfit_f <- muhaz2(Surv(surv_mm/12, status == 1) ~ 1, data = melanoma)
## scale hazard by 1000
plot(hazfit_f, xlab="Years since diagnosis",col="blue",lty="solid", haz.scale=1000, xlim=c(0,10))

##(g) Estimate incidence rate ratios as a function of follow-up 
poisson_g <- glm( status == 1 ~ fu + offset( log(pt) ),
                                     family = poisson,
                                     data = melanoma_10y_spl )
summary(poisson_g)
eform(poisson_g)

##(h) Poisson regression adjusting for time since diagnosis
poisson_h <- glm( status == 1 ~ fu + year8594 + offset( log(pt) ),
                          family = poisson,
                          data = melanoma_10y_spl )
summary(poisson_h)
eform(poisson_h)

# Add interaction
poisson_h2 <- glm( status == 1 ~ fu*year8594 + offset( log(pt) ), 
                   family=poisson, 
                   data=melanoma_10y_spl )

summary(poisson_h2)
eform(poisson_h2)

##(i)-i. Poisson regression adjusting for age, calendar period and sex
poisson_i <- glm( status == 1 ~ fu + year8594 + sex + agegrp + offset(log(pt)), 
                   family=poisson, 
                   data=melanoma_10y_spl)
summary(poisson_i)
eform(poisson_i)

##(i)-ii.  Test if the effect of agegrp is significant using a likelihood ratio test
drop1(poisson_i, ~agegrp, test="Chisq")

##(i)-iii. Wald test of the overall effect of age and interpret the results
linearHypothesis(poisson_i,c("agegrp1 = 0","agegrp2 = 0","agegrp3 = 0"))
# Explanation: agegrp is highly significant in the model.

##(j) Is the effect of sex modified by calendar period
poisson_j <- glm( status == 1 ~ fu + agegrp + year8594*sex + offset( log(pt) ), 
                  family=poisson, data=melanoma_10y_spl )
summary(poisson_j)
eform(poisson_j)
# Explanation: The interaction term (year8594Diagnosed 85-94:sex2) is not 
# statistically significant indicating the effect of sex is modified by period.

##(k)-i. Calculate effect of sex for year8594==2
hr_k <- exp(coef(poisson_j))
hr_k["sex2"]
# Explanation: The effect of sex for patients diagnosed 1975–84 is 0.6031338.

hr_k["sex2"]*hr_k["year8594Diagnosed 85-94:sex2"]
# Explanation: The effect of sex for patients diagnosed 1985–94 is 0.5691922.
