setwd("~/Documents/survival/")
library(survival)
library(survminer)
library(tidyverse)

data <- read.table("eBLvssBL_survival.txt", sep="\t", header=TRUE)

fit <- survfit(Surv(OS_TTE_censored, OS_Event) ~cohort, data=data)
print(fit)

summary(fit)


ggsurvplot(fit=fit,
           data=data,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_bw(),
           palette = c("#e51b18", "#377eb8"))
