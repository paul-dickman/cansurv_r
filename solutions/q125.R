## Exercise 125
## Created: 2021-04-09 Enoch Chen 
## Edited:  2021-04-09 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q22.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr)

## Read data
brv <- read_dta("brv.dta")

#' variables are coded as follows
#' couple: id for couple
#' dob: Date of birth
#' doe: Date of entry into study
#' dox: Date of exit from study
#' dosp:  Date of bereavement
#' fail:  0 = alive, 1=died  
#' group: 1, 2, 3
#' disab: 0, 1, 2, 3
#' health: Perceived health status: 0, 1, 2
#' sex: 1=M, 2=F

## Check the loaded data
head(brv)

##(a) list selected variables
brv %>% select(id, sex, doe, dosp, dox, fail, couple) %>%
        filter(couple== 3 | couple== 4 | couple== 19 | couple == 7)  %>%
        arrange(couple, id)

##(b) Crude mortality rate for each sex
scale <- 365.24
brv <- mutate(brv, age_entry = as.numeric(doe - dob) / scale, # Calc age at entry
                   att_age = as.numeric(dox - dob) / scale,   # Calc attained age
                   t_at_risk = att_age - age_entry)           # Calc time at risk

## Crude rates per 1000 person-years
survRate(Surv(age_entry/1000, att_age/1000, fail) ~ sex, data=brv)

## Poisson regression
poisson_b <- glm( fail ~ sex + offset( log( t_at_risk) ), family=poisson, data=brv)
summary(poisson_b)
eform(poisson_b)

## Calculate mean entry age by sex
select(brv, sex, age_entry) %>%
       group_by(sex) %>%
       summarise(mean_age_entry = mean(age_entry))

##(c) Breaking records into pre and post bereavement
## Creating times relative to spouse death (year=0)
brv <- mutate(brv,
               id=NULL,
               y_before_sp_dth =  as.numeric(doe -dosp) / scale,
               y_after_sp_dth = as.numeric(dox - dosp) / scale)

## Splitting at spouse death (year=0)
brv_split <- survSplit(brv, cut = 0, end="y_after_sp_dth", start="y_before_sp_dth", id="id",event="fail")

## Calculating risk times
brv_split <- mutate(brv_split,
                   t_sp_at_risk =   y_after_sp_dth - y_before_sp_dth,
                   brv = ifelse(y_after_sp_dth > 0, 1, 0))

## Take a look at the five first couples
brv_split %>% select(couple, id, sex, doe, dosp, dox, fail, y_before_sp_dth, y_after_sp_dth, t_sp_at_risk)  %>% 
              filter(couple== 3 | couple== 4 | couple== 19 | couple == 7)  %>%
              arrange(couple, id)

##(d) Poisson regression
poisson_d <- glm( fail ~ brv + offset( log(t_sp_at_risk) ), family=poisson, data=brv_split)
summary(poisson_d)
eform(poisson_d)

##(e) Poisson regression stratified on sex
## sex = 1, Male
poisson_e1 <- glm( fail ~ brv + offset( log(t_sp_at_risk) ), family=poisson, data=filter(brv_split, sex==1))
summary(poisson_e1)
eform(poisson_e1)

## sex = 2, Female
poisson_e2 <- glm( fail ~ brv + offset( log(t_sp_at_risk) ), family=poisson, data=filter(brv_split, sex==2))
summary(poisson_e2)
eform(poisson_e2)

## Poisson regression with effect of bereavement by sex 
brv_split <- mutate(brv_split,
                    sex = as.factor(sex),
                    brv = as.factor(brv))

poisson_e3 <- glm( fail ~ sex + brv:sex + offset( log(t_sp_at_risk) ), family=poisson, data=brv_split)
summary(poisson_e3)
eform(poisson_e3)

##(f) Controlling for age.
## Look at mean and lowest/highest value for _t0 and _t
brv_split <- brv_split %>%
             mutate(age_sp_dth =  as.numeric(dosp - dob) / scale, # Age at spouse death
             t0 = age_sp_dth + y_before_sp_dth,      # Age at start of timeband
             t  = age_sp_dth + y_after_sp_dth)       # Age at end of timeband

summarise(brv_split, mean(t0), min(t0), max(t0))
summarise(brv_split, mean(t), min(t), max(t))

## Split follow up time (=age) to be able to control for age bands in Poisson model
age_cat <- seq(70,100,5) 

brv_split <- survSplit(brv_split, cut=age_cat, start="t0", end="t", event="fail", zero = 0)

## Creating age band category
brv_split <- mutate(brv_split,
                    age = cut(t, age_cat),
                    t_at_risk = t - t0 )  

## Crude rates
survRate(Surv(t_at_risk, fail) ~ age, data=brv_split)

## Poisson model controlling for age
poisson_f <- glm( fail ~ brv + age + offset( log(t_at_risk) ), family=poisson, data=brv_split)
summary(poisson_f)
eform(poisson_f)

## Poisson regression controlling for age and sex 
poisson_f2 <- glm( fail ~ brv + age + sex + offset( log(t_at_risk) ), family=poisson, data=brv_split)
summary(poisson_f2)
eform(poisson_f2)

##(g) Poisson regression estimating the effect of brv for each sex, controlling for age 
poisson_g <- glm( fail ~ age + sex + brv:sex + offset( log(t_at_risk) ), family=poisson, data=brv_split)
summary(poisson_g)
eform(poisson_g)

##(h)
#' No syntax needed

##(i) Cox regression estimate the effect of bereavement adjusted for age
summary(coxph(Surv(t0, t, fail) ~ brv,
              data = brv_split))

summary(coxph(Surv(t0, t, fail) ~ brv + sex ,
              data = brv_split))

##(j) Cox regression estimating the effect of brv for each sex, controlling for age
summary(coxph(Surv(t0, t, fail) ~ sex + sex:brv,
              data = brv_split))