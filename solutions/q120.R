## Exercise 120
## Created: 2021-03-26 Enoch Chen 
## Edited:  2021-03-26 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q9.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr)   # manipulate data
library(splines) # natural splines
library(car)     # linearHypothesis

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

# Restrict follow-up to 10 years
melanoma_10y <- mutate(melanoma,
                       status= ifelse(surv_mm <= 120 & status == 1, 1, 0),
                       surv_mm = ifelse(surv_mm<=120, surv_mm, 120) # censored everyone after 120 months
)

##(a) Cox regression
coxfit_a <- coxph(Surv(surv_mm, status == 1) ~ year8594,
               data = melanoma_10y)
summary(coxfit_a)

##(b) Log-rank test
survdiff(Surv(surv_mm, status == 1) ~ year8594, 
         data = melanoma_10y)

##(c) Cox regression including sex and age
coxfit_c <- coxph(Surv(surv_mm, status == 1) ~ year8594 + sex + agegrp,
                           data = melanoma_10y)
summary(coxfit_c)

##(c)-i. For patients of the same sex diagnosed in the same calendar period, 
#        agegrp 2 at diagnosis have an estimated 86% higher risk of death due to melanoma 
#        than those agegrp 0 at diagnosis.

##(c)-ii. The parameter estimate for period changes from 0.78 to 0.72 
#         when age and sex are added to the model. 

##(c)-iii. Wald test
linearHypothesis(coxfit_c,c("agegrp1 = 0","agegrp2 = 0","agegrp3 = 0"))

# agegrp is highly significant in the model.

##(d) Test if the effect of agegrp is significant
coxfit_d <- coxph(Surv(surv_mm, status == 1) ~ year8594 + sex,
                  data = melanoma_10y)
summary(coxfit_d)

## LR test
anova(coxfit_c, coxfit_d)

# agegrp is highly significant in the model.
# The Wald test is an approximation to the LR test and we would expect the two to be similar.

##(e)-i.
# coxfit_c <- coxph(Surv(surv_mm, status == 1) ~ year8594 + sex + agegrp,
#                   data = melanoma_10y)

# When splitting the Poisson regression model 
# we split time since diagnosis into annual intervals and explicitly estimated the effect of this
# factor in the model. The Cox model does not estimate the effect of `time' but the other
# estimates are adjusted for `time'.

# Poisson regression
# Split the time scale
melanoma_10y_spl <- survSplit(melanoma_10y, cut=12*(0:10), end="surv_mm", start="start",
                              event = "status" ) # event must be a charater variable

## Calculate follow-up time and recode start time as a factor
melanoma_10y_spl <- mutate(melanoma_10y_spl,
                           pt = surv_mm - start,
                           fu = as.factor(start) )

# Poisson regression adjusted by year8594, sex, agegrp
# Same as q111(i)-i.
poisson_e <- glm( status == 1 ~ fu + year8594 + sex + agegrp + offset(log(pt)), 
                  family=poisson, 
                  data=melanoma_10y_spl)
summary(poisson_e)
eform(poisson_e)

##(e)-ii. Compare with Cox regression
summary(coxfit_c)

##(e)-iii. 
#' Yes, both models assume ‘proportional hazards’. 
#' The proportional hazards assumption implies that the risk ratios 
#' for sex, period, and age are constant across all levels of follow-up time. 
#' In other words, the assumption is that there is no effect modification by follow-up time.

##(f) Please refer to the results from (e)

##(g) By splitting at each failure time we can estimate a Poisson regression model 
#     that is identical to the Cox model
fit_g <- survfit(Surv(surv_mm, event=status == 1) ~ 1,
                  data = melanoma_10y )

## Have a look at the fitted object
str(fit_g, 1)
head(fit_g$time)

## Split on time since diagnosis (1-year intervals)
melanoma_10y_spl2 <- survSplit(melanoma_10y, cut=fit_g$time, 
                               end="surv_mm", start="start",
                               event="status")

## Calculate follow-up time and recode start time as a factor
melanoma_10y_spl2 <- mutate(melanoma_10y_spl2,
                           pt = surv_mm - start,
                           fu = as.factor(start) )

## Collapse
melanoma_10y_spl2 <- melanoma_10y_spl2 %>%
                     group_by(fu,year8594,sex,agegrp) %>%
                     summarise(pt=sum(pt), death_cancer=sum(status == 1)) %>%
                     data.frame

## Poisson regression
poisson_g <- glm( death_cancer ~ fu + year8594 + sex + agegrp + offset( log(pt) ),
                  family = poisson,
                  data = melanoma_10y_spl2 )
summary(poisson_g)
eform(poisson_g)

##(h) Split the data finely (e.g., 3-month intervals) and model the effect of time 
# using a restricted cubic spline with Poisson regression. 

## Split the time scale by 3 months 
cuts.splines <- seq(0, max(fit_g$time), by=3)
mid.splines <- cuts.splines + 1.5
melanoma_10y_spl3 <- survSplit(Surv(surv_mm, status == 1)~., data=melanoma_10y, 
                               cut=cuts.splines, start="tstart", end="tstop") %>%
                     mutate(cut = cut(tstop, cuts.splines),
                            mid = mid.splines[unclass(cut)]) %>%
                     group_by(mid, year8594, sex, agegrp) %>%
                     summarise(pt = sum(tstop-tstart), status = sum(event))

poisson_h <- glm( status ~ ns(mid, df=3) + year8594 + sex + agegrp + offset( log(pt) ),
                  family = poisson,
                  data = melanoma_10y_spl3 )
summary(poisson_h)
eform(poisson_h)

## Quick approach: use the effects package
library(effects)
plot(Effect("mid", poisson_h))

## Slow approach do it step by step
## Utility function to draw a confidence interval
polygon.ci <- function(time, interval, col="lightgrey") 
               polygon(c(time,rev(time)), 
                       c(interval[,1],rev(interval[,2])), 
                       col=col, border=col)

## Define exposures
times <- seq(0,max(cuts.splines),length=1000)
delta <- diff(times)[1]
newdata <- data.frame(mid=times+delta/2, year8594="Diagnosed 85-94",
                      sex="1", agegrp="2", # Male, agegrp 45-59
                      pt=1)

## Predict rates and 95% CI 
## The example in Stata plots log hazard, here we plot hazard
pred <- predict(poisson_h, newdata=newdata, se.fit=TRUE)
predrate <- exp(pred$fit)
ci <- with(pred, exp(cbind(fit-1.96*se.fit, fit+1.96*se.fit)))

##Plot
matplot(newdata$mid, ci, type="n", xlab="Time since diagnosis (months)",
        ylab="Rate", main="Males aged 45-59 years diagnosed 1985-94")
polygon.ci(newdata$mid, ci) 
lines(newdata$mid, predrate)
