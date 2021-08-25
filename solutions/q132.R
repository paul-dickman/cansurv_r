## Exercise 132
## Created: 2021-08-05 Enoch Chen 
## Edited:  2021-08-10 Enoch Chen
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
years <- levels(melanoma$year8594)
agegrps <- levels(melanoma$agegrp)

newdata1 <- data.frame(female=0, year8594=years[1], agegrp2=0, agegrp3=0, agegrp4=0)
newdata2 <- data.frame(female=0, year8594=years[1], agegrp2=1, agegrp3=0, agegrp4=0)
newdata3 <- data.frame(female=0, year8594=years[1], agegrp2=0, agegrp3=1, agegrp4=0)
newdata4 <- data.frame(female=0, year8594=years[1], agegrp2=0, agegrp3=0, agegrp4=1)

plot(fpm_b, newdata=newdata1, 
     xlab="Time since diagnosis (years)",
     ylab="Cause specific mortality rate (per pys)",
     type="haz", ci=FALSE, ylim=c(0,0.2))

lines(fpm_b, newdata=newdata2, type="haz", lty=2)
lines(fpm_b, newdata=newdata3, type="haz", lty=3)
lines(fpm_b, newdata=newdata4, type="haz", lty=4)
legend("topright", legend=paste0(agegrps), lty=1:4)
 
##(c) Time-dependent effects for age group
fpm_c <- stpm2(Surv(t, event) ~ female + year8594 + agegrp2 + agegrp3 + agegrp4, 
               tvc=list(agegrp2= 2 , agegrp3= 2 , agegrp4 = 2), data = melanoma, df=4)
summary(fpm_c)
eform(fpm_c)

#Likelihood-ratio test
anova(fpm_b, fpm_c)

##(d) predict the hazard for each age group
plot(fpm_b, newdata=newdata1, 
     xlab="Time since diagnosis (years)",
     ylab="Cause specific mortality rate (per pys)",
     type="haz", ci=FALSE, ylim=c(0,0.2), line.col = "red")
lines(fpm_c, newdata=newdata1, type="haz", lty=2, col = "red")

lines(fpm_b, newdata=newdata2, type="haz", lty=1, col = "blue")
lines(fpm_c, newdata=newdata2, type="haz", lty=2, col = "blue")

lines(fpm_b, newdata=newdata3, type="haz", lty=1, col = "purple")
lines(fpm_c, newdata=newdata3, type="haz", lty=2, col = "purple")

lines(fpm_b, newdata=newdata4, type="haz", lty=1, col = "green")
lines(fpm_c, newdata=newdata4, type="haz", lty=2, col = "green")

legend("topright", legend=paste0(agegrps), lty=1, col=c("red","blue","purple","green"))

##(e) time-dependent hazard ratios
plot(fpm_c, newdata=newdata1, 
     xlab="Time since diagnosis (years)",
     ylab="Log hazard ratio",
     log="y",
     type="hr", var="agegrp2", ci=FALSE, line.col = "blue", ylim=c(0, 50))

lines(fpm_c, newdata=newdata1, type="hr", var="agegrp3" , col = "red")
lines(fpm_c, newdata=newdata1, type="hr", var="agegrp4" , col = "green")

plot(fpm_c, newdata=newdata1,
     type="hr", log="y",
     var="agegrp4", 
     xlab="Time since diagnosis (years)",
     ylab="Log hazard ratio"
     )

##(f) Difference in hazard rates
plot(fpm_c,newdata=newdata1,
     type="hdiff", var="agegrp4",
     xlab="Time since diagnosis (years)")

##(g) Difference in survival functions
newdata_g1 <- data.frame(female=1, year8594=years[2], agegrp2=0, agegrp3=0, agegrp4=0)
newdata_g2 <- data.frame(female=1, year8594=years[2], agegrp2=0, agegrp3=0, agegrp4=1)

# Plot survival function of agegrp1 and agegrp2 separately
plot(fpm_c,newdata=newdata_g1,
     type="surv", 
     xlab="Time since diagnosis (years)", ci = FALSE, ylim = c(0.8,1) ,
     line.col = "blue")
lines(fpm_c,newdata=newdata_g2,
      type="surv", col = "red")

legend("topright", legend=c("Age <45", "Age 75+"), lty=1, col=c("blue","red"))

# Plot survival difference
plot(fpm_c, newdata = newdata_g1,
     type="sdiff", var="agegrp4",
     xlab="Time since diagnosis (years)",
     ylab="Difference in survival functions",
     ylim = c(-0.3, 0))

## (h) varying df for time-dependent effects
lapply(1:3, function(i){
  hr4_df_i <- stpm2(Surv(t, event) ~ female + year8594 + agegrp2 + agegrp3 + agegrp4, 
                  tvc=list(agegrp2= i , agegrp3= i , agegrp4 = i), data = melanoma, df=4)
  data.frame(
    i,
    AIC=AIC(hr4_df_i),
    BIC=BIC(hr4_df_i))
    
})

# tvc = 1
hr4_df1 <- stpm2(Surv(t, event) ~ female + year8594 + agegrp2 + agegrp3 + agegrp4, 
                    tvc=list(agegrp2= 1 , agegrp3= 1 , agegrp4 = 1), data = melanoma, df=4)

# tvc = 2
hr4_df2 <- stpm2(Surv(t, event) ~ female + year8594 + agegrp2 + agegrp3 + agegrp4, 
                 tvc=list(agegrp2= 2 , agegrp3= 2 , agegrp4 = 2), data = melanoma, df=4)

# tvc = 3
hr4_df3 <- stpm2(Surv(t, event) ~ female + year8594 + agegrp2 + agegrp3 + agegrp4, 
                 tvc=list(agegrp2= 3 , agegrp3= 3 , agegrp4 = 3), data = melanoma, df=4)

# Plot 3 lines in one figure
plot(hr4_df1, newdata=newdata1, 
     xlab="Time since diagnosis (years)",
     ylab="Log hazard ratio",
     log = "y",
     type="hr", var="agegrp4", ci=FALSE, line.col = "red", ylim=c(0.5, 500))

lines(hr4_df2, newdata=newdata1, type="hr", var="agegrp4" , col = "blue")
lines(hr4_df3, newdata=newdata1, type="hr", var="agegrp4" , col = "green")

legend("topright", legend=c("1 df", "2 df", "3 df"), lty=1, col=c("red","blue","green"))

##(i) two time-dependent effects
fpm_i <- stpm2(Surv(t, event) ~ female + agegrp2 + agegrp3 + agegrp4, 
                 tvc=list(agegrp2= 3, agegrp3= 3, agegrp4 = 3, female = 3), data = melanoma, df=4)
summary(fpm_i)
eform(fpm_i)

newdata_i1 <- data.frame(female=0, agegrp2=0, agegrp3=0, agegrp4=0)
newdata_i2 <- data.frame(female=0, agegrp2=0, agegrp3=0, agegrp4=1)

plot(fpm_i, newdata=newdata_i1, 
     xlab="Time since diagnosis (years)",
     type="hr", 
     var="female",
     ci=TRUE, ci.col = rgb(red = 1, green = 0, blue = 0, alpha = 0.1), line.col = "red")

lines(fpm_i, newdata=newdata_i2, 
      type="hr", 
      var="female",
      ci=TRUE, ci.col = rgb(red = 0, green = 0, blue = 1, alpha = 0.1) , col = "blue")

##(j) 
##'skip
##'Please refer to the demonstration in Stata



