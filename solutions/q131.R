## Exercise 131
## Created: 2021-06-21 Enoch Chen 
## Edited:  2021-06-21 Enoch Chen
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

##(b) Weibull model (using stpm2)z
fpm_b <- stpm2(Surv(t, event) ~ 1, data = melanoma, df=1)
summary(fpm_b)
eform(fpm_b)

s <- predict(fpm_b, grid=FALSE, full=TRUE, se.fit=FALSE, type="surv")

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

