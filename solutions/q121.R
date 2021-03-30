## Exercise 121
## Created: 2021-03-30 Enoch Chen 
## Edited:  2021-03-30 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q10.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr)   # manipulate data
library(ggplot2)
#library(splines) # natural splines
#library(car)     # linearHypothesis

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
                       surv_mm = ifelse(surv_mm<=120, surv_mm, 120), # censored everyone after 120 months
                       surv_yy = surv_mm/12,
                       )

##(a) Hazard function by year8594
haz_a <- muhaz2(Surv(surv_yy, status == 1) ~ year8594, data = melanoma_10y)
haz_a_df <- as.data.frame(haz_a)

## Plot hazard
plot(haz_a, haz.scale=1000,
     xlab="Time since diagnosis (years)", 
     ylab="Hazard per 1000 person-years")

## Alternative: ggplot2
ggplot(haz_a_df, aes(x=x, y=y*1000, colour= strata)) + geom_line() +
  xlab("Time since diagnosis (years)") +
  ylab("Hazard per 1000 person-years")

##(b) Hazard function on log scale 
plot(haz_a, log="y")

##(c) Log cumulative hazard function
survfit_c <- survfit(Surv(surv_yy, status == 1) ~ year8594, data=melanoma_10y)
plot(survfit_c, col=1:2, fun=function(S) - log(-log(S)), log="x",
     xlab="log(time)", ylab="-log(H)")

##Alternative: biostat3::survPHplot
biostat3::survPHplot(Surv(surv_yy, status == 1) ~ year8594, data=melanoma_10y)

##(d) Cox regression
cox_d <- coxph(Surv(surv_yy, status == 1) ~ year8594, data=melanoma_10y)
summary(cox_d)

##(e) Cox regression adjusted for sex, calendar period and age 
cox_e <- coxph(Surv(surv_yy, status == 1) ~ sex + year8594 + agegrp, data=melanoma_10y)
summary(cox_e)

##Plot of the scaled Schoenfeld residuals for calendar period 1985–94
cox_e.phtest <- cox.zph(cox_e, transform="identity") 
plot(cox_e.phtest[2], resid=TRUE, se=TRUE, main="Test of PH assumption: Schoenfeld residuals", ylim=c(-4,4))

##(f) Redo (a)-(e) for agegrp
##Hazard function by agegrp
haz_f <- muhaz2(Surv(surv_yy, status == 1) ~ agegrp, data = melanoma_10y)
haz_f_df <- as.data.frame(haz_f)

## Plot hazard
plot(haz_f, haz.scale=1000,
     xlab="Time since diagnosis (years)", 
     ylab="Hazard per 1000 person-years")

## Alternative: ggplot2
ggplot(haz_f_df, aes(x=x, y=y*1000, colour= strata)) + geom_line() +
  xlab("Time since diagnosis (years)") +
  ylab("Hazard per 1000 person-years")

##Hazard function on log scale 
plot(haz_f, log="y")

## Log cumulative hazard function
survfit_f <- survfit(Surv(surv_yy, status == 1) ~ agegrp, data=melanoma_10y)
plot(survfit_f, col=1:4, fun=function(S) - log(-log(S)), log="x",
     xlab="log(time)", ylab="-log(H)")
legend("topright",legend=levels(melanoma_10y$agegrp),col=1:4,lty=1)

##Alternative: biostat3::survPHplot
biostat3::survPHplot(Surv(surv_yy, status == 1) ~ agegrp, data=melanoma_10y)

## Cox regression
cox_f <- coxph(Surv(surv_yy, status == 1) ~ agegrp, data=melanoma_10y)
summary(cox_f)

## Cox regression adjusted for sex, calendar period and age 
cox_f <- coxph(Surv(surv_yy, status == 1) ~ sex + year8594 + agegrp, data=melanoma_10y)
summary(cox_f)

##Plot of the scaled Schoenfeld residuals for calendar period 1985–94
cox_f.phtest <- cox.zph(cox_f, transform="identity") 
plot(cox_f.phtest[1], resid=TRUE, se=TRUE, main="Test of PH assumption: Schoenfeld residuals", ylim=c(-4,4))
plot(cox_f.phtest[2], resid=TRUE, se=TRUE, main="Test of PH assumption: Schoenfeld residuals", ylim=c(-4,4))
plot(cox_f.phtest[3], resid=TRUE, se=TRUE, main="Test of PH assumption: Schoenfeld residuals", ylim=c(-4,4))

## (g) Formally test the PH assumption
print(cox_f.phtest)

## (h) Split the time scale by follow-up time
melanoma_10y_spl <- survSplit(melanoma_10y, cut=c(2), end="surv_yy", start="start",
                              event="status", episode="fu") %>%
                              mutate(fu = as.factor(fu))

##Alternative 1: Cox regression using time-dependent option 
cox_h <- coxph(Surv(start, surv_yy, status == 1) ~ sex + year8594 + agegrp*fu, data=melanoma_10y_spl)
summary(cox_h)

##Alternative2: Split data
##Tabulate ageband including risk_time
melanoma_10y_spl %>% select(id, start, surv_yy) %>% filter(id<=3) %>% arrange(id, surv_yy)

##(i) Cox regression with interaction