## Exercise 131
## Created: 2021-06-21 Enoch Chen 
## Edited:  2021-07-02 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q28.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr)   # manipulate data
library(rstpm2)  # nsx
library(ggplot2) # ggplot

## Read data
melanoma <- read_dta("melanoma.dta")%>% 
            filter(stage == 1) %>%
            mutate(year8594 = ifelse(year8594 == 1, "Diagnosed 85-94", "Diagnosed 75-84"),
                   female = ifelse(sex == 2, 1, 0), # Create a female variable
                   event  = ifelse((status == 1) & surv_mm<= 120.5, 1, 0), # Redefine status
                   surv_mm = ifelse(surv_mm<= 120.5, surv_mm, 120.5),     # censored everyone after 120 months
                   t = surv_mm/12,                                     # surv_mm in years
                   agegrp1 = (agegrp== 1)+0, # used by time-dependent effect
                   agegrp2 = (agegrp== 2)+0, # used by time-dependent effect
                   agegrp3 = (agegrp== 3)+0, # used by time-dependent effect
                   agegrp4 = (agegrp== 4)+0) # used by time-dependent effect
                       
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

##(a) Kaplan-Meier curve
km <- survfit(Surv(t, event == 1) ~ 1, data = melanoma) # make Kaplan-Meier estimates
summary(km)      

# Plot Kaplan-Meier curve
plot(km,             
     ylab="S(t)",
     xlab="analysis time",
     main = "Kaplan−Meier survival estimate")

##(b) Weibull model (using stpm2)
fpm_b <- stpm2(Surv(t, event) ~ 1, data = melanoma, df=1)
summary(fpm_b)
eform(fpm_b)

# levels() can only be used in factor variables
time <- levels(as.factor(melanoma$t))

s <- predict(fpm_b, grid=TRUE,full=TRUE,se.fit=TRUE,
             newdata = data.frame(t = time),
             type="surv")

ggplot(s, aes(x=t,y=Estimate, ymin=0.6, ymax=1)) +
       xlab("analysis time") +
       ylab("Survival") + 
       geom_line() 

## Overlay KM curve with stpm2-Weibull
plot(km,             
     ylab="S(t)",
     xlab="analysis time",
     main = "Kaplan−Meier survival estimate")
lines(s$t, s$Estimate, lwd = 1)

##(c) Plot hazard function
h <- predict(fpm_b, grid=TRUE,full=TRUE,se.fit=TRUE,
             newdata = data.frame(t = time),
             type="haz")

# To plot smoothed hazards, we use the muhaz package
plot(muhaz2(Surv(t, event) ~ 1, data = melanoma),
     xlim = c(0,10),
     ylim = c(0,0.05))

# Add stpm2-fitted hazard
lines(h$t, h$Estimate, lwd = 1, col = "red")
legend("bottomright", legend = c("Smoothed hazard function", "stpm2"), col = c("black", "red"), lty=c("solid", "solid"))

##(d)  try stpm2 with 4 df
fpm_d <- stpm2(Surv(t, event) ~ 1, data = melanoma, df=4)
summary(fpm_d)
eform(fpm_d)

s <- predict(fpm_d, grid=TRUE,full=TRUE,se.fit=TRUE,
             newdata = data.frame(t = time),
             type="surv")

ggplot(s, aes(x=t,y=Estimate, ymin=0.6, ymax=1)) +
        xlab("analysis time") +
        ylab("Survival") + 
        geom_line() 

## Survival function
plot(km,             
     ylab="S(t)",
     xlab="analysis time",
     main = "Kaplan−Meier survival estimate")
lines(s$t, s$Estimate, lwd = 1,  col = "red")

legend("bottomright", legend = c("Survival function", "stpm2"), col = c("black", "red"), lty=c("solid", "solid"))

## Hazard function
h <- predict(fpm_d, grid=TRUE,full=TRUE,se.fit=TRUE,
             newdata = data.frame(t = time),
             type="haz")

# To plot smoothed hazards, we use the muhaz package
plot(muhaz2(Surv(t, event) ~ 1, data = melanoma),
     xlim = c(0,10),
     ylim = c(0,0.05))

# Add stpm2-fitted hazard
lines(h$t, h$Estimate, lwd = 1, col = "red")
legend("bottomright", legend = c("Smoothed hazard function", "stpm2"), col = c("black", "red"), lty=c("solid", "solid"))

##(e) Fit a Cox Model
cox_e <- coxph(Surv(t, event) ~ year8594, data = melanoma)
summary(cox_e)

##(f) Equivalent flexible parametric model 
fpm_f <- stpm2(Surv(t, event) ~ year8594, data = melanoma, df=4)
summary(fpm_f)
eform(fpm_f)

##(e) Predicted survival and hazard functions by period of diagnosis
years <- levels(as.factor(melanoma$year8594))

## Survival function
s <- predict(fpm_f,newdata=data.frame(year8594=years),
             grid=TRUE,full=TRUE,se.fit=TRUE,
             type="surv")

head(s)

ggplot(s, aes(x=t, y=Estimate, fill = year8594, ymin=lower, ymax=upper)) +
        xlab("Time since diagnosis (years)") +
        ylab("Survival") +
        geom_ribbon(alpha=0.6) +
        geom_line()

## Hazard function
plot(fpm_f,newdata=data.frame(year8594=years[1]), type="haz", ci=FALSE,
     xlab="Time since diagnosis (years)", ylab="Cause-specific mortality rate (per py)")
lines(fpm_f,newdata=data.frame(year8594=years[2]), type="haz", col=2)
legend("topright", legend=years, col=1:2, lty=1, bty="n")

##(h) hazard on log scale
plot(fpm_f,newdata=data.frame(year8594=years[1]), type="haz", ci=FALSE,
     xlab="Time since diagnosis (years)", ylab="Log scale of cause-specific mortality rate (per py)", log="y")
lines(fpm_f,newdata=data.frame(year8594=years[2]), type="haz", col=2)
legend("topright", legend=years, col=1:2, lty=1, bty="n")

##(i) sensitivity to knots
Rbind <- function(object) do.call(rbind, object)
out <- lapply(1:6, function(i) {
        fitaic <- stpm2(Surv(t, event) ~ year8594, data=melanoma, df=i)
        data.frame(
                i,
                AIC=AIC(fitaic),
                BIC=BIC(fitaic),
                beta=as.numeric(coef(fitaic)[2]),
                se=coef(summary(fitaic))[2,2])
})
out %>% Rbind

##(j) Compare hazard and survival
## Baseline survival
fitaic0 <- stpm2(Surv(t, event) ~ year8594, data=melanoma, df=6)
plot(fitaic0,newdata=data.frame(year8594=years[1]), lty=6, ci=FALSE,
     xlab="Time since diagnosis (years)")
for (i in 1:5 ) {
        fitaic <- stpm2(Surv(t, event) ~ year8594, data=melanoma, df=i)
        lines(fitaic,newdata=data.frame(year8594=years[1]), lty=i)
}
legend("topright", legend=paste0("df=",1:6), lty=1:6)

## Baseline hazard
fitaic1 <- stpm2(Surv(t, event) ~ year8594, data=melanoma, df=6)
plot(fitaic1,newdata=data.frame(year8594=years[1]), lty=6, type="haz",
     ci=FALSE, xlab="Time since diagnosis (years)", ylab="Hazard")
for (i in 1:5 ) {
        fitaic <- stpm2(Surv(t, event) ~ year8594, data=melanoma, df=i)
        lines(fitaic,type="haz",newdata=data.frame(year8594=years[1]), lty=i)
}
legend("topright", legend=paste0("df=",1:6), lty=1:6)

##(k) knot locations
##'run the following code to fit 10 models with 5df (6 knots) where
##'the 4 internal knots are selected at random centiles of the 
##'distribution of event times.

## Specily the spline placement
knots <- quantile(melanoma$t, probs=c(0, .2, .4, .6, .8))

melanoma$spline <- nsx(melanoma$t, knots = knots)

fpm_k <- stpm2(Surv(t, event) ~ year8594 + spline, data=melanoma)
fpm_k <- stpm2(Surv(t, event) ~ year8594 , data=melanoma, df =5 )

summary(fpm_k)

##'As there are many ties in this data we add a small random number to the survival times
##'(otherwise we risk having knots in the same location)
melanoma$t2 <- melanoma$t + rnorm(nrow(melanoma), mean = 0, sd = 0.001)

fpm_k <- stpm2(Surv(t2, event) ~ year8594, data=melanoma, df=c(22, 33))

plot(fpm_k,newdata=data.frame(year8594=years[1]), lty=6, ci=FALSE,
     xlab="Time since diagnosis (years)")
for (i in 1:5 ) {
        fpm_k <- stpm2(Surv(t2, event) ~ year8594, data=melanoma, df=i)
        lines(fitaic,newdata=data.frame(year8594=years[1]), lty=i)
}
legend("topright", legend=paste0("df=",1:6), lty=1:6)
