## Exercise 112
## Created: 2021-03-22 Enoch Chen 
## Edited:  2021-03-22 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q8.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr) # manipulate data
library(bshazard)  # for nonparametric smoothing of hazard function

## Read data
diet <- read_dta("diet.dta")

#' hieng (energy intake) is coded as follows
#' 0 [low]
#' 1 [high]

diet <- mutate(diet, hieng = as.factor(ifelse(hieng == 1, "high", "low")),
                     y1k = y / 1000)

## Check the loaded data
head(diet)

##(a) Incidence rates by attained age
scale <- 365.24

# Incidence of hieng == "low"
plot(bshazard(Surv((doe - dob)/scale, (dox - dob)/scale, chd) ~ 1, 
              data=subset(diet, hieng == "low")),
     ylim=c(0,0.03), conf.int=FALSE, xlab="Attained age (years)")

# Incidence of hieng == "high"
lines(bshazard(Surv((doe - dob)/scale, (dox - dob)/scale, chd) ~ 1, 
               data=subset(diet, hieng == "high")),
      col="red", conf.int=FALSE)

legend("topleft", legend=c('hieng=="low"','hieng=="high"'), col=1:2, lty=1, bty="n")

# Incidence rates by time since diagnosis
plot(muhaz2(Surv((dox - doe)/scale, chd) ~ hieng, data=diet), lty=1,
     xlab="Time since study entry (years)", ylim=c(0,0.025),
     legend=FALSE)

legend("topleft", legend=c('hieng=="low"','hieng=="high"'), col=1:2, lty=1, bty="n")


##(b) Modelling the rate, no adjustment for timescale
# Make the reference group hieng = "low"
diet$hieng <- relevel(diet$hieng, ref = "low")

poisson_b <- glm( chd ~ hieng + offset( log( y1k ) ), 
                  family=poisson, data=diet)

summary(poisson_b)
eform(poisson_b)

##(c) Adjustment for confounders job and bmi 
# Create BMI variable
diet$bmi <- diet$weight/((diet$height/100)^2)

#' job is coded as follows
#' 1 [driver]
#' 2 [conductor]
#' 3 [bank]

poisson_c <- glm( chd ~ hieng + job + bmi + offset( log( y1k ) ), 
                  family=poisson, data=diet)
summary(poisson_c)
eform(poisson_c)

##(d) 
