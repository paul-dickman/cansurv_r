## Excercise 101
## Created: 2021-03-02 Enoch Chen 
## Edited:  2021-03-02 Enoch Chen 
###############################################################################
## Load the biostat3 library
library(biostat3)

## Life table
print(lifetab2(Surv(floor(surv_yy), status == "Dead: cancer")~1, colon_sample, breaks=0:10), digits=2)
print(lifetab2(Surv(floor(surv_mm), status == "Dead: cancer")~1, colon_sample, breaks=c(12*(0:10))), digits=2)

## Make Kaplan-Meier estimates
## Dead: cancer is Status == 1
mfit <- survfit(Surv(surv_mm, status == "Dead: cancer") ~ 1, data = colon_sample)
summary(mfit)          

## Kaplan-Meier plot 
plot(mfit,                                            
ylab="S(t)",
xlab="Time since diagnosis in months",
main = "Kaplan−Meier estimates of cause−specific survival")

## Kaplan-Meier plot with risk table

## Load survminer
library(survminer)

ggsurvplot(
  mfit,                      # survfit object.
  risk.table = TRUE,         # show risk table.
  conf.int = TRUE,           # show confidence intervals for
  break.time.by = 20       , # numeric value controlling time axis breaks
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE  # show bars instead of names in text annotations
  # in legend of risk table
)