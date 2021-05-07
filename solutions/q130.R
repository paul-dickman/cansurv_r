## Exercise 130
## Created: 2021-04-16 Enoch Chen 
## Edited:  2021-05-07 Enoch Chen
###############################################################################
## Load the packages
library(biostat3)
library(haven)
library(dplyr)   # manipulate data
library(splines) # natural splines
library(ggplot2)

## Read data
## All stages
melanoma <- read_dta("melanoma.dta")%>% 
            mutate(year8594 = ifelse(year8594 == 1, "Diagnosed 85-94", "Diagnosed 75-84"),
                   agegrp = as.factor(agegrp),
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
                       status  = ifelse((status == 1 | status == 2) & surv_mm<= 120, 1, 0), # Redefine status
                       surv_mm = ifelse(surv_mm<= 120, surv_mm, 120),     # censored everyone after 120 months
                       t = surv_mm/12                                     # surv_mm in years
                       )

##(a)-i. Split the data with narrow (1 month) time intervals
# Split the time scale
melanoma_10y_spl <- survSplit(melanoma_10y, cut=1/12*(0:120), end="t", start="start", event = "status" ) # event must be a character variable

## Calculate follow-up time and recode start time as a factor
melanoma_10y_spl <- mutate(melanoma_10y_spl,
                           pt = t - start,
                           fu = as.factor(start) )
## Collapse
melanoma_10y_spl2 <- melanoma_10y_spl %>%
                     group_by(fu, female, year8594, agegrp) %>%
                     summarise(pt=sum(pt), d = sum(status == 1)) 

## Generate interval 
melanoma_10y_spl2 <- melanoma_10y_spl2  %>%
                     group_by(fu) %>%
                     mutate(interval = group_indices()) ## equivalent to Stata's egen interval=group()

##(a)-ii.Fit a model with a parameter for each interval
## glm -1 is suppress constant term, equivalent to Stata's nocons
poisson_a <- glm( d ~ as.factor(interval) -1 + offset( log(pt)), 
                  family = poisson,
                  data = melanoma_10y_spl2 )
summary(poisson_a)
eform(poisson_a)

##  predict the baseline (one parameter for each interval)
melanoma_10y_spl3 <- melanoma_10y_spl2 %>% mutate (pt0 = pt, ## keep the original pt
                                                   pt = 1 )  ## to have log(pt) = 0 as no offset
melanoma_10y_spl3$haz_grp <- predict(poisson_a, newdata = melanoma_10y_spl3,
                                     type = "response")
## per 1000 person-years
melanoma_10y_spl3$haz_grp_1k <- melanoma_10y_spl3$haz_grp * 1000

## Plot
ggplot(data = melanoma_10y_spl3, aes(x=interval/12, y = haz_grp_1k)) +
    geom_point()+
    scale_y_continuous(name="Baseline hazard (1000 pys)", breaks = c(5,10,20,50,100,150), limits = c(0, 150))+
    scale_x_continuous(name="Years from diagnosis", breaks = c(2, 4, 6, 8, 10), limits = c(0, 10)) +
    ggtitle("Localised skin melanoma. Plot of the estimated baseline hazard function for the piecewise
model.")
  

##(b) linear splines (1 knot at knots at 1.5 years)
melanoma_10y_spl4 <- melanoma_10y_spl3 %>% 
                     mutate(midtime = 0.0416667 + 0.0833333*(interval-1),
                            lin_s1 = midtime,
                            lin_int2 = ifelse(midtime>1.5, 1, 0),
                            lin_s2 = (midtime - 1.5)*lin_int2)

## Fit two separate linear regression lines (4 parameters)
poisson_b <- glm( d ~ lin_s1 + lin_int2 + lin_s2 + offset( log(pt0)), 
                  family = poisson,
                  data = melanoma_10y_spl4 )

summary(poisson_b)
eform(poisson_b)

##  predict the baseline (one parameter for each interval)
melanoma_10y_spl4 <- melanoma_10y_spl4 %>% mutate ( pt1 = pt0, ## keep the original pt
                                                    pt0 = 1 )  ## to have log(pt) = 0 as no offset
melanoma_10y_spl4$haz_lin1 <- predict(poisson_b, newdata = melanoma_10y_spl4,
                                     type = "response")
## per 1000 person-years
melanoma_10y_spl4$haz_lin1_1k <- melanoma_10y_spl4$haz_lin1 * 1000

## Plot
ggplot(data = melanoma_10y_spl4, aes(x=midtime, y = haz_grp_1k)) +
  geom_point()+
  geom_line(data = melanoma_10y_spl4%>%filter(midtime<=1.5), aes(x=midtime, y = haz_lin1_1k), color = "red") +
  geom_line(data = melanoma_10y_spl4%>%filter(midtime>1.5), aes(x=midtime, y = haz_lin1_1k), color = "red") +
  geom_vline(xintercept = 1.5, linetype="dashed") +
  scale_y_continuous(name="Baseline hazard (1000 pys)", breaks = c(5,10,20,50,100,150), limits = c(0, 150))+
  scale_x_continuous(name="Years from diagnosis", breaks = c(2, 4, 6, 8, 10), limits = c(0, 10)) +
  ggtitle("Localised skin melanoma. Plot of the estimated baseline hazard function for the piecewise
model and linear spline model.")

## The gradient up to 1.5 years is: _b[lin_s1]
summary(poisson_b)$coefficients[2]

## The gradient after 1.5 years is: _b[lin_s1] + _b[lin_s2]
summary(poisson_b)$coefficients[2] + summary(poisson_b)$coefficients[4]


##(c) Force the functions to join at the knot (3 parameters)		
poisson_c <- glm( d ~ lin_s1 + lin_s2 + offset( log(pt1)), 
                  family = poisson,
                  data = melanoma_10y_spl4 )

summary(poisson_c)
eform(poisson_c)

##  predict the baseline (one parameter for each interval)
melanoma_10y_spl4 <- melanoma_10y_spl4 %>% mutate ( pt2 = pt1, ## keep the original pt
                                                    pt1 = 1 )  ## to have log(pt) = 0 as no offset
melanoma_10y_spl4$haz_lin2 <- predict(poisson_c, newdata = melanoma_10y_spl4,
                                      type = "response")
## per 1000 person-years
melanoma_10y_spl4$haz_lin2_1k <- melanoma_10y_spl4$haz_lin2 * 1000

## Plot
ggplot(data = melanoma_10y_spl4, aes(x=midtime, y = haz_grp_1k)) +
  geom_point()+
  geom_line(data = melanoma_10y_spl4, aes(x=midtime, y = haz_lin2_1k), color = "red") +
  geom_vline(xintercept = 1.5, linetype="dashed") +
  scale_y_continuous(name="Baseline hazard (1000 pys)", breaks = c(5,10,20,50,100,150), limits = c(0, 150))+
  scale_x_continuous(name="Years from diagnosis", breaks = c(2, 4, 6, 8, 10), limits = c(0, 10)) +
  ggtitle("Localised skin melanoma. Plot of the estimated baseline hazard function for the piecewise
model and linear spline model.")

## The gradient up to 1.5 years is: _b[lin_s1]
summary(poisson_c)$coefficients[2]

## The gradient after 1.5 years is: _b[lin_s1] + _b[lin_s2]
summary(poisson_c)$coefficients[2] + summary(poisson_c)$coefficients[3]

## (d) Now use cubic polynomials with 1 knot at 2 years
melanoma_10y_spl5 <- melanoma_10y_spl4 %>% 
                     mutate(cubic_s1 = midtime,
                            cubic_s2 = midtime^2,
                            cubic_s3 = midtime^3,
                            cubic_int = ifelse(midtime>2, 1, 0), 
                            cubic_lin = (midtime - 2)*cubic_int,
                            cubic_quad = ((midtime - 2)^2)*cubic_int,
                            cubic_s4 = ((midtime - 2)^3)*cubic_int)

## Regression
poisson_d <- glm( d ~ cubic_s1 + cubic_s2 + cubic_s3 + cubic_int + cubic_lin+ cubic_quad + cubic_s4+ offset( log(pt2)), 
                  family = poisson,
                  data = melanoma_10y_spl5 )

summary(poisson_d)
eform(poisson_d)

##  predict the baseline (one parameter for each interval)
melanoma_10y_spl5 <- melanoma_10y_spl5 %>% mutate ( pt3 = pt2, ## keep the original pt
                                                    pt2 = 1 )  ## to have log(pt) = 0 as no offset
melanoma_10y_spl5$haz_cubic1 <- predict(poisson_d, newdata = melanoma_10y_spl5,
                                      type = "response")
## per 1000 person-years
melanoma_10y_spl5$haz_cubic1_1k <- melanoma_10y_spl5$haz_cubic1 * 1000

## Plot
ggplot(data = melanoma_10y_spl5, aes(x=midtime, y = haz_grp_1k)) +
  geom_point()+
  geom_line(data = melanoma_10y_spl5%>%filter(midtime<=2), aes(x=midtime, y = haz_cubic1_1k), color = "red") +
  geom_line(data = melanoma_10y_spl5%>%filter(midtime>2), aes(x=midtime, y = haz_cubic1_1k), color = "red") +
  geom_vline(xintercept = 2, linetype="dashed") +
  scale_y_continuous(name="Baseline hazard (1000 pys)", breaks = c(5,10,20,50,100,150), limits = c(0, 150))+
  scale_x_continuous(name="Years from diagnosis", breaks = c(2, 4, 6, 8, 10), limits = c(0, 10)) +
  ggtitle("Localised skin melanoma. Plot of the estimated baseline hazard function for the piecewise
model and cubic spline model.")

## (e) constrain to join at knots (drop separate intercept)	
## Regression
poisson_e <- glm( d ~ cubic_s1 + cubic_s2 + cubic_s3 + cubic_lin+ cubic_quad + cubic_s4+ offset( log(pt3)), 
                  family = poisson,
                  data = melanoma_10y_spl5 )

summary(poisson_e)
eform(poisson_e)

##  predict the baseline (one parameter for each interval)
melanoma_10y_spl5 <- melanoma_10y_spl5 %>% mutate ( pt4 = pt3, ## keep the original pt
                                                    pt3 = 1 )  ## to have log(pt) = 0 as no offset
melanoma_10y_spl5$haz_cubic2 <- predict(poisson_e, newdata = melanoma_10y_spl5,
                                        type = "response")
## per 1000 person-years
melanoma_10y_spl5$haz_cubic2_1k <- melanoma_10y_spl5$haz_cubic2 * 1000

## Plot
ggplot(data = melanoma_10y_spl5, aes(x=midtime, y = haz_grp_1k)) +
  geom_point()+
  geom_line(data = melanoma_10y_spl5, aes(x=midtime, y = haz_cubic2_1k), color = "red") +
  geom_vline(xintercept = 2, linetype="dashed") +
  scale_y_continuous(name="Baseline hazard (1000 pys)", breaks = c(5,10,20,50,100,150), limits = c(0, 150))+
  scale_x_continuous(name="Years from diagnosis", breaks = c(2, 4, 6, 8, 10), limits = c(0, 10)) +
  ggtitle("Localised skin melanoma. Plot of the estimated baseline hazard function for the piecewise
model and cubic spline model.")


## (f) continuous 1st derivative (drop second linear term)
poisson_f <- glm( d ~ cubic_s1 + cubic_s2 + cubic_s3 + cubic_quad + cubic_s4 + offset( log(pt4)), 
                  family = poisson,
                  data = melanoma_10y_spl5 )

summary(poisson_f)
eform(poisson_f)

##  predict the baseline (one parameter for each interval)
melanoma_10y_spl5 <- melanoma_10y_spl5 %>% mutate ( pt5 = pt4, ## keep the original pt
                                                    pt4 = 1 )  ## to have log(pt) = 0 as no offset
melanoma_10y_spl5$haz_cubic3 <- predict(poisson_f, newdata = melanoma_10y_spl5,
                                        type = "response")
## per 1000 person-years
melanoma_10y_spl5$haz_cubic3_1k <- melanoma_10y_spl5$haz_cubic3 * 1000

## Plot
ggplot(data = melanoma_10y_spl5, aes(x=midtime, y = haz_grp_1k)) +
  geom_point()+
  geom_line(data = melanoma_10y_spl5, aes(x=midtime, y = haz_cubic3_1k), color = "red") +
  geom_vline(xintercept = 2, linetype="dashed") +
  scale_y_continuous(name="Baseline hazard (1000 pys)", breaks = c(5,10,20,50,100,150), limits = c(0, 150))+
  scale_x_continuous(name="Years from diagnosis", breaks = c(2, 4, 6, 8, 10), limits = c(0, 10)) +
  ggtitle("Localised skin melanoma. Plot of the estimated baseline hazard function for the piecewise
model and cubic spline model with continuous first derivatives.")

## (g) continuous 2nd derivative (drop second quadratic term)
poisson_g <- glm( d ~ cubic_s1 + cubic_s2 + cubic_s3 + cubic_s4 + offset( log(pt5)), 
                  family = poisson,
                  data = melanoma_10y_spl5 )

summary(poisson_g)
eform(poisson_g)

##  predict the baseline (one parameter for each interval)
melanoma_10y_spl5 <- melanoma_10y_spl5 %>% mutate ( pt6 = pt5, ## keep the original pt
                                                    pt5 = 1 )  ## to have log(pt) = 0 as no offset
melanoma_10y_spl5$haz_cubic4 <- predict(poisson_g, newdata = melanoma_10y_spl5,
                                        type = "response")
## per 1000 person-years
melanoma_10y_spl5$haz_cubic4_1k <- melanoma_10y_spl5$haz_cubic4 * 1000

## Plot
ggplot(data = melanoma_10y_spl5, aes(x=midtime, y = haz_grp_1k)) +
  geom_point()+
  geom_line(data = melanoma_10y_spl5, aes(x=midtime, y = haz_cubic4_1k), color = "red") +
  geom_vline(xintercept = 2, linetype="dashed") +
  scale_y_continuous(name="Baseline hazard (1000 pys)", breaks = c(5,10,20,50,100,150), limits = c(0, 150))+
  scale_x_continuous(name="Years from diagnosis", breaks = c(2, 4, 6, 8, 10), limits = c(0, 10)) +
  ggtitle("Localised skin melanoma. Plot of the estimated baseline hazard function for the piecewise
model and cubic spline model with continuous first and second derivatives.")

##restricted cubic splines
library(Epi) ## use natural splines from Epi package

##(h) generate splines with 5 knots (4 df) ?????
fit <- lm(d ~ ns(midtime, df = 4), data = melanoma_10y_spl5)
summary(fit)
attr(terms(fit), "predvars")

cbind(melanoma_10y_spl5, Ns(melanoma_10y_spl5$midtime, df=4))

