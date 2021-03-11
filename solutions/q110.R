## Exercise 110
## Created: 2021-03-11 Enoch Chen 
## Edited:  2021-03-11 Enoch Chen
## Reference: Biostatistics III in R. https://biostat3.net/download/R/solutions/q6.html
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr) # manipulate data

## Read data
diet <- read_dta("diet.dta")

## Check the loaded data
head(diet)

##(a) Tabulate CHD incidence rates per 1000 person-years for each category of hieng.
# Generate a new time variable y1k
diet$y1k <- diet$y/1000

## Estimated rate (per 1000 person-years)
diet %>%
     group_by(hieng) %>%
     summarise(Event = sum(chd), Time = sum(y1k), Rate = Event/Time,     
               CI_low = poisson.test(Event,Time)$conf.int[1],
               CI_high = poisson.test(Event,Time)$conf.int[2])

## Caculate rate ratio
7.07/13.6

## The same output can be generated as follows 
## diet.irr <- survRate(Surv(y1k,chd) ~ hieng, data=diet)
## with(diet.irr, poisson.test(event,tstop))           #Denominator: hieng == 1
## with(diet.irr, poisson.test(rev(event),rev(tstop))) #Denominator: hieng == 0

##(b) Poisson regression for incidence rate ratios (high energy to low)
## glm (outcome ~ exposure + offset (time))
poisson_b <- glm( chd ~ hieng + offset( log( y1k ) ), family=poisson, data=diet)
summary(poisson_b)

## Be aware of taking the exponential of the estimates
eform(poisson_b)

##(c) Histogram 
hist_c <- hist(diet$energy, breaks=25, probability=TRUE, xlab="Total energy (kcals per day)")
curve(dnorm(x, mean=mean(diet$energy), sd=sd(diet$energy)), col = "red", add=TRUE)

## Summarise
quantile(diet$energy, probs=c(0.01,0.05,0.1,0.25,0.5,0.75,0.90,0.95,0.99))

##(d) Generate new variable eng3
diet$eng3 <- cut(diet$energy, breaks=c(1500,2500,3000,4500),
                 labels=c("low","medium","high"), 
                 right = FALSE)

## Make table
cbind(Freq=table(diet$eng3),
      Prop=table(diet$eng3)/nrow(diet))

##(e) Estimate and plot the rates for difference levels of eng3
diet_irr_e <- survRate(Surv(y/1000,chd) ~ eng3, data=diet)
print(diet_irr_e)

## IRR with CIs
with(diet_irr_e, rate[eng3=="medium"] / rate[eng3=="low"])
        with(diet_irr_e[c(2,1),], { # compare second row with first row
        poisson.test(event, tstop)})
        
with(diet_irr_e, rate[eng3=="high"] / rate[eng3=="low"])
        with(diet_irr_e[c(3,1),], { # compare third row with first row
        poisson.test(event, tstop)})
        
##(f) Create your own indicator variables for the three levels of eng3 
## Make dummy variables for eng3
diet <- mutate(diet, 
               X1 = as.numeric(eng3 == "low"),
               X2 = as.numeric(eng3 == "medium"),
               X3 = as.numeric(eng3 == "high"))

## Check the first 5 rows of each group
head(filter(diet, eng3=="low") %>% 
     select(c("eng3", "X1", "X2", "X3")), 5)

head(filter(diet, eng3=="medium") %>% 
     select(c("eng3", "X1", "X2", "X3")), 5)

head(filter(diet, eng3=="medium") %>%  
     select(c("eng3", "X1", "X2", "X3")), 5)

##(h) Poisson regression to compare second and third to first
poisson_h <- glm( chd ~ X2 + X3 + offset( log( y1k ) ), family=poisson, data=diet )
summary(poisson_h)
eform(poisson_h)

# Explanation: The estimated IRRs with 95% CIs
# X2 to X1: 0.6452 (0.3389 - 1.2284)
# X3 to X1: 0.2886 (0.1235 - 0.6744)

##(i) Repeat (h) but using level 2 as reference level
poisson_i <- glm( chd ~ X1 + X3 + offset( log( y1k ) ), family=poisson, data=diet )
summary(poisson_i)
eform(poisson_i)

# Explanation: The estimated IRRs with 95% CIs
# X1 to X2: 1.5498 (0.8140 - 2.9506)
# X3 to X2: 0.4473 (0.1991 - 1.0047)

##(j) Repeat (i) but using automatically created indicator
## In R, if a variable is a factor, 
## R would automatically make dummy variables of it.
class(diet$eng3) # eng3 is factor

poisson_j <- glm( chd ~ eng3 + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson_j )
eform( poisson_j )

# Estimates and explantion are identical to (h)

##(k) Calculate the total number of events during follow-up (person-time at risk)
## That is, per person-year (rather than per 1000 person-years)
summarise(diet, failures = sum(chd), 
                total_person_time = sum(y),
                rate = sum(chd) / sum(y))

# The estimated incidence rate is 0.00999 events per person-year
