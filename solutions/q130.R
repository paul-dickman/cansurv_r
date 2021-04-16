## Exercise 130
## Created: 2021-04-16 Enoch Chen 
## Edited:  2021-04-16 Enoch Chen
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr)   # manipulate data
library(splines) # natural splines


## Read data
## All stages
melanoma <- read_dta("melanoma.dta")%>% 
            mutate(year8594 = ifelse(year8594 == 1, "Diagnosed 85-94", "Diagnosed 75-84"),
                   agegrp = as.factor(agegrp),
                   sex = as.factor(sex),
                   female = ifelse(sex == 2, 1, 0),
                   stage = as.factor(stage),
                   subsite = as.factor(subsite))

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
#' 
#' subsite is coded as follows
#' 1  [Head and Neck]
#' 2  [Trunk]
#' 3  [Limbs]
#' 4  [Multiple and NOS]
#' 
#' stage is coded as follows
#' 0  [Unknown]
#' 1  [Localised]
#' 2  [Regional] 
#' 3  [Distant]

# Restrict follow-up to 10 years and create new variable for all-cause death and cause-specific death
melanoma_10y <- mutate(melanoma,
                       status  = ifelse(status == 1 | status == 2, 1, 0),# Redefine status
                       surv_mm = ifelse(surv_mm<=120, surv_mm, 120),     # censored everyone after 120 months
                       t = surv_mm/12                                    # surv_mm in years
                       )

## Split the data with narrow (1 month) time intervals
# Split the time scale
melanoma_10y_spl <- survSplit(melanoma_10y, cut=1/12*(0:120), end="t", start="start", event = "status" ) # event must be a character variable

## Calculate follow-up time and recode start time as a factor
melanoma_10y_spl <- mutate(melanoma_10y_spl,
                           pt = t - start,
                           fu = as.factor(start) )
## Collapse
melanoma_10y_spl2 <- melanoma_10y_spl %>%
                     group_by(fu, female, year8594, agegrp) %>%
                     summarise(pt=sum(pt), d = sum(status == 1)) %>%
                     data.frame %>%
                     group_by(fu) %>%
                     mutate(interval = group_indices()) ## equivalent to Stata's egen interval=group()

## TODO: Why there were more people died at the last interval??
## The model looks right but output is different from Stata.
##(a) Fit a model with a parameter for each interval
poisson_a <- glm( d ~ as.factor(interval) + offset( log(pt) ),
                  family = poisson,
                  data = melanoma_10y_spl2 )
summary(poisson_a)
eform(poisson_a)
