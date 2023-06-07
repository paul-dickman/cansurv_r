## Exercise 231
## Created: 2023-06-07 Enoch Chen
## Notes: only finished q231 (a) and (b)
###############################################################################
# Load packages
library(haven)
library(dplyr)
library(rstpm2)
library(ggplot2)

##(a)
# Load colon dta
colon <- read_dta("colon.dta") 
popmort <- read_dta("popmort.dta") 

# Prepare colon dataset
colon <- colon %>% 
         mutate(status= ifelse(surv_mm <= 60.5 & status %in% c(1, 2), 1, 0),
                surv_mm = ifelse(surv_mm <= 60.5, surv_mm, 60.5),
                st_years = surv_mm/12,
                att_age = pmin(floor(age + st_years), 99),
                att_year = floor(yydx + st_years),           
                sex = as.numeric(sex))

# Merge with popmort file
colon <- colon %>% 
         left_join(popmort,
                   by = c("sex" = "sex",
                          "att_age" = "_age",
                          "att_year" = "_year")) %>%
         mutate(female = ifelse(sex == 2, 1, 0)) %>%
         filter(age<=90)

##(b)
# Flexible parametric relative survival model
fit <- stpm2(Surv(st_years, status == 1) ~ age, 
             data = colon, df = 5,
             bhazard(rate))
summary(fit)
eform(fit)

# Predict the hazard functions for the specified age values 50, 60, 70
# At time 0 to 5 with intervals 0.5
exh <- predict(fit, type = "hazard", 
               newdata = expand.grid(age = c(50, 60, 70), st_years=seq(0, 5, by = 0.05)),
               full=TRUE, se.fit=TRUE)
exh$Estimate_1000 <- exh$Estimate*1000
  
ggplot(exh, aes(x=st_years, y=Estimate_1000, color=factor(age), group=age)) +
       labs(x = "Time", y = "Excess hazard (per 1000 person-years)", 
            title = "Hazard Functions for Different Ages") +
       geom_line()

                 