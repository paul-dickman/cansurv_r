## Exercise 112
## Created: 2021-03-22 Enoch Chen 
## Edited:  2021-03-23 Enoch Chen
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
# Create BMI variable & make job into factor variable
diet <- diet%>%mutate(bmi = weight/(height/100)^2,
                      job = as.factor(job))

#' job is coded as follows
#' 1 [driver]
#' 2 [conductor]
#' 3 [bank]

poisson_c <- glm( chd ~ hieng + job + bmi + offset( log( y1k ) ), 
                  family=poisson, data=diet)
summary(poisson_c)
eform(poisson_c)

##(d)  Modelling the rate, adjusting for timescale age
# Split ageband
# Split time at 30,50,60 and 72 with time scale age at entry to attained age
scale <- 365.24
ageband <- c(30,50,60,72)
diet.spl <- survSplit(Surv((doe - dob)/scale, (dox - dob)/scale, chd) ~ ., 
                          data=diet, cut=ageband, start="t0",end="t")

# Tabulate ageband
diet.spl %>% select(id, t0, t, y) %>% 
         slice(c(1:10)) %>% 
         arrange(id, t0)

# Generate a risk time variable
diet.spl$risktime <- diet.spl$t - diet.spl$t0

# Tabulate ageband including risk_time
diet.spl %>% select(id, t0, t, y, risktime) %>% 
             slice(c(1:10)) %>% 
             arrange(id, t0)

# Tabulate ageband chd
diet.spl$agespan <- cut(diet.spl$t, ageband)
xtabs(~agespan+chd,diet.spl)

# Poisson regression adjusted for ageband
# Affirm that the reference group is hieng = low
diet.spl <- within(diet.spl, hieng <- relevel(hieng, ref = "low"))
poisson_d <- glm( chd ~ hieng + agespan + offset( log( risktime) ),
                  family=poisson,
                  data=diet.spl)

summary(poisson_d)
eform(poisson_d)

# Adjustment for confounders job, bmi
poisson_d2 <- glm( chd ~ hieng + agespan + job + bmi + offset( log( risktime) ),
                  family=poisson,
                  data=diet.spl)

summary(poisson_d2)
eform(poisson_d2)

##(e) Modelling the rate, adjusting for timescale time-in-study
fuband <- c(0, 5, 10, 15, 22)
diet.spl.t_entry <- survSplit(Surv((dox-doe)/365.24, chd) ~ ., 
                              data=diet, cut=fuband, end="t", start="t0", 
                              event="chd")

##Tabulate fuband
diet.spl.t_entry %>% select(id, t0, t, y) %>% 
                     slice(c(1:10)) %>% 
                     arrange(id, t0)

#Tabulate fuband including risktime
diet.spl.t_entry$risktime <- diet.spl.t_entry$t - diet.spl.t_entry$t0

diet.spl.t_entry %>% select(id, t0, t, y, risktime) %>% 
                     slice(c(1:10)) %>% 
                     arrange(id, t0)

# Tabulate fuband chd
diet.spl.t_entry$fuspan <- cut(diet.spl.t_entry$t, fuband)
xtabs(~fuspan+chd,diet.spl.t_entry)

# Poisson regression adjusted for fuband
# Affirm that the reference group is hieng = low
diet.spl.t_entry <- within(diet.spl.t_entry, hieng <- relevel(hieng, ref = "low"))
poisson_e <- glm( chd ~ hieng + fuspan + offset( log( risktime) ),
                  family=poisson,
                  data=diet.spl.t_entry)

summary(poisson_e)
eform(poisson_e)

# Adjustment for confounders job, bmi
poisson_e2 <- glm( chd ~ hieng + fuspan + job + bmi + offset( log( risktime) ),
                   family=poisson,
                   data=diet.spl.t_entry)

summary(poisson_e2)
eform(poisson_e2)
