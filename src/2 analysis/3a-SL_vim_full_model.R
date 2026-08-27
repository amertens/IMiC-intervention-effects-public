# 3a-SL_vim_full_model.R
# =============================================================================
# SL arm-classification VIM pipeline, split 1 of 3 (was 3-SL_vim.R).
#
# Fits the FULL cross-validated SuperLearner classifier (all milk-component
# biomarkers as covariates, predicting treatment arm), study x visit x
# contrast. This is the fastest of the three stages and the one shared
# dependency the other two (3b, 3c) build on -- `full_AUC_res` gets combined
# into both the individual-lab and leave-one-out result tables.
#
# Output: data/models/SL_vim_full_fit.RDS  (list(full_AUC_res, full_SL_fits))
#
# Run from the repo root:
#   Rscript "src/2 analysis/3a-SL_vim_full_model.R"
# =============================================================================

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

library(r2weight)
library(quadprog)
library(tidyverse)
library(ck37r)
library(data.table)

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d <- readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

#expand the dataset with contrast variable to group_by()
table(d$studyid, d$arm)
d1 <- d %>% filter(studyid=="ELICIT", arm %in% c("Control","Az.")) %>% mutate(contrast="Az.")
d2 <- d %>% filter(studyid=="ELICIT", arm %in% c("Control","Nico")) %>% mutate(contrast="Nico")
d3 <- d %>% filter(studyid=="ELICIT", arm %in% c("Control","Nico+Az.")) %>% mutate(contrast="Nico+Az.")
d4 <- d %>% filter(studyid=="MISAME-3", arm %in% c("Control","BEP/BEP")) %>% mutate(contrast="BEP/BEP")
d5 <- d %>% filter(studyid=="MISAME-3", arm %in% c("Control","BEP/IFA")) %>% mutate(contrast="BEP/IFA")
d6 <- d %>% filter(studyid=="MISAME-3", arm %in% c("Control","IFA/BEP")) %>% mutate(contrast="IFA/BEP")
d7 <- d %>% filter(studyid=="VITAL-Lactation", arm %in% c("Control","BEP+ExBf")) %>% mutate(contrast="BEP+ExBf")
d8 <- d %>% filter(studyid=="VITAL-Lactation", arm %in% c("Control","BEP+ExBf+AZT")) %>% mutate(contrast="BEP+ExBf+AZT")
d <- bind_rows(d1,d2,d3,d4,d5,d6,d7,d8)

Xvars = as.vector(unlist(all_milk_components))

table(d$arm)
table(d$contrast)

#### full fit model (all biomarkers) ------------------------------------------
set.seed(12345)
full_SL_fits = fit_SuperLearner(d,outcome="arm", covars=as.vector(unlist(all_milk_components)), slmod=sl, CV=T, family="binomial")
full_AUC_res = extract_AUC(full_SL_fits)
full_AUC_res

saveRDS(list(full_AUC_res=full_AUC_res,
             full_SL_fits=full_SL_fits),
        file=paste0(here::here(),"/data/models/SL_vim_full_fit.RDS"))

cat("wrote data/models/SL_vim_full_fit.RDS\n")
