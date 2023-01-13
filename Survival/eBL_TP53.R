setwd("~/Documents/survival/")
library(survival)
library(survminer)
library(tidyverse)

# Load data
data <- read.table("eBL_TP53_surv.txt", sep="\t", header=TRUE)

# Mutations
fit <- survfit(Surv(OS_time_censored, OS_event_censored) ~TP53_mut, data=data)
print(fit)

summary(fit)
ggsurvplot(fit=fit,
           data=data,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Overall survival analysis of eBL patients stratified by TP53 mutation status")


fit <- survfit(Surv(RR_time_censored, RR_event_censored) ~TP53_mut, data=data)
print(fit)

summary(fit)


ggsurvplot(fit=fit,
           data=data,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Risk of relapse analysis of eBL patients stratified by TP53 mutation status")

# All data OS
fit <- survfit(Surv(OS_time_censored, OS_event_censored) ~TP53_any, data=data)
print(fit)

summary(fit)
ggsurvplot(fit=fit,
           data=data,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Overall survival analysis of eBL patients stratified by the presence of any TP53 abnormality")


# All data PFS
fit <- survfit(Surv(RR_time_censored, RR_event_censored) ~TP53_any, data=data)
print(fit)

summary(fit)
ggsurvplot(fit=fit,
           data=data,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Risk of relapse analysis of eBL patients stratified by the presence of any TP53 abnormality")

