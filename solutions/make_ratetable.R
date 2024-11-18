## Make ratetable
## Created: 2023-06-05 Joshua Entrop 
## Edited:  2023-06-06 Enoch Chen: minor changes
## Updated: 2024-11-18 Paul Dickman
###############################################################################
## Load the packages
library(haven)
library(tidyverse)

## Read the (Stata format) popmort file
popmort <- read_dta("popmort.dta")

## Prepare popmort file for making ratetable
popmort_wide <- pivot_wider(popmort,
                            id_cols = c("sex", "_age"),
                            names_from  = "_year",
                            values_from = "prob")

popmort_males <- popmort_wide %>% 
                 filter(sex == 1)           %>% 
                 select(-c("sex", "_age"))  %>% 
                 as.matrix()

popmort_females <- popmort_wide %>% 
                   filter(sex == 2)          %>% 
                   select(-c("sex", "_age")) %>% 
                   as.matrix()

ratetable <- transrate(popmort_males,
                       popmort_females,
                       yearlim = c(1951, 2000),
                       int.length = 1)

