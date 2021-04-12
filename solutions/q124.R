## Exercise 124
## Created: 2021-04-07 Enoch Chen 
## Edited:  2021-04-07 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q13.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)

## Read data
diet <- read_dta("diet.dta")

## Check the loaded data
head(diet)

##(a) Poisson regression
poisson_a <- glm( chd ~ hieng + offset( log( y ) ), family=poisson, data=diet)
summary(poisson_a)
eform(poisson_a)   # IRR

## Time in study as timescale
## Cox regression
cox_a <- coxph(Surv(y, chd) ~ hieng, data=diet)
summary(cox_a)

##(b) Attained age as timescale 
## Cox regression
scale <- 365.24
cox_b <- coxph(Surv((doe-dob)/scale, (dox-dob)/scale, chd) ~ hieng, data=diet)
summary(cox_b)