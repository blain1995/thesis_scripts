setwd("~/Documents/survival/")
library(survival)
library(survminer)
library(tidyverse)

# Load data
data <- read.table("sBL_TP53_AUG.txt", sep="\t", header=TRUE)
sanger <- subset(data, Alex_paper == 1)

# Overall Survival
fit <- survfit(Surv(OS_TTE_censored_3yr, OS_Event_3yr) ~TP53_mut, data=sanger)
print(fit)

summary(fit)

ggsurvplot(fit=fit,
           data=sanger,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Overall survival analysis of sBL patients stratified by TP53 mutation status")


fit <- survfit(Surv(OS_TTE_censored_3yr, OS_Event_3yr) ~TP53_CN, data=sanger)
print(fit)

summary(fit)

ggsurvplot(fit=fit,
           data=sanger,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Overall survival analysis of sBL patients stratified by TP53 copy number status")

fit <- survfit(Surv(OS_TTE_censored_3yr, OS_Event_3yr) ~TP53_LOH, data=sanger)
print(fit)

summary(fit)

ggsurvplot(fit=fit,
           data=sanger,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Overall survival analysis of sBL patients stratified by TP53 LOH status")


fit <- survfit(Surv(OS_TTE_censored_3yr, OS_Event_3yr) ~TP53_abn, data=sanger)
print(fit)

summary(fit)

ggsurvplot(fit=fit,
           data=sanger,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Overall survival analysis of sBL patients stratified by the presence any TP53 abnormality")

# Progression-free Survival
fit <- survfit(Surv(PFS_TTE_3yr, PFS_Event_3yr) ~TP53_mut, data=sanger)
print(fit)

summary(fit)

ggsurvplot(fit=fit,
           data=sanger,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Progression-free analysis of sBL patients stratified by TP53 mutation status")


fit <- survfit(Surv(PFS_TTE_3yr, PFS_Event_3yr) ~TP53_CN, data=sanger)
print(fit)

summary(fit)

ggsurvplot(fit=fit,
           data=sanger,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Progression-free analysis of sBL patients stratified by TP53 copy number status")

fit <- survfit(Surv(PFS_TTE_3yr, PFS_Event_3yr) ~TP53_LOH, data=sanger)
print(fit)

summary(fit)

ggsurvplot(fit=fit,
           data=sanger,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Progression-free analysis of sBL patients stratified by TP53 LOH status")


fit <- survfit(Surv(PFS_TTE_3yr, PFS_Event_3yr) ~TP53_abn, data=sanger)
print(fit)

summary(fit)

ggsurvplot(fit=fit,
           data=sanger,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Progression-free analysis of sBL patients stratified by the presence any TP53 abnormality")


fit <- survfit(Surv(PFS_TTE_3yr, PFS_Event_3yr) ~TP53_mut, data=sanger)
print(fit)

summary(fit)

ggsurvplot(fit=fit,
           data=sanger,
           pval=TRUE,
           conf.int=FALSE,
           risk.table = TRUE,
           risk.table.col = "strata",
           ggtheme= theme_classic(),
           palette = c("#377eb8", "#e51b18"),
           legend="none", 
           title="Progression-free analysis of sBL patients stratified by TP53 mutation status")
