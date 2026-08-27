
# SUPERSEDED 2026-08-25: split into three independently-runnable stages so a
# crash partway through (this script died mid-run repeatedly on background
# Windows profile switches) doesn't lose completed work and each stage can be
# rerun on its own:
#   3a-SL_vim_full_model.R      -- full-model fit (fast; the shared dependency)
#   3b-SL_vim_individual_lab.R  -- single-group fits (feeds manuscript Fig 1)
#   3c-SL_vim_leave_one_out.R   -- leave-one-group-out fits (not on the Fig 1 path)
# This file is kept for reference; do not run it going forward.

#------------------------------------------------------------------------------
# SL cvAUC for variable importance in predicting arm
#------------------------------------------------------------------------------


rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

library(r2weight)
library(quadprog)
library(tidyverse)
library(ck37r)
library(data.table)

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d <- readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))


#need to add microbiome and untargeted metabolomics

#update to test all arms

# d <- d %>% filter( arm %in% c("Control","BEP/BEP","Nico+Az.","BEP+ExBf+AZT")) %>% mutate(arm=ifelse(arm %in% c("BEP/BEP","Nico+Az.","BEP+ExBf+AZT"),1,0))
# table(d$arm)

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

#run full fit model 
set.seed(12345)
full_SL_fits = fit_SuperLearner(d,outcome="arm", covars=as.vector(unlist(all_milk_components)), slmod=sl, CV=T, family="binomial")
full_AUC_res = extract_AUC(full_SL_fits)
full_AUC_res

saveRDS(list(full_AUC_res=full_AUC_res, 
             full_SL_fits=full_SL_fits),
        file=paste0(here::here(),"/data/models/SL_vim_full_fit.RDS"))



#loop over sets of biomarkers
SL_fits_list <- list()
for(i in 1:length(all_milk_components)){
  #print(i)
  SL_fits=AUC_res=NULL
  set.seed(12345)
  SL_fits = fit_SuperLearner(d,outcome="arm", covars=Xvars[!(Xvars %in% all_milk_components[[i]])], slmod=sl, CV=T, family="binomial")
  AUC_res = extract_AUC(SL_fits)
  AUC_res
  
  SL_fits_list[[i]] <- list(SL_fits=SL_fits, AUC_res=AUC_res)
}
names(SL_fits_list) <- names(all_milk_components)
saveRDS(SL_fits_list,file=paste0(here::here(),"/data/models/SL_vim_fits.RDS"))

#save combined AUC tables for plotting
full_AUC_res <- full_AUC_res %>% mutate(group='all')
vim_bio <- rbindlist(lapply(SL_fits_list, function(x) x[[2]]), idcol='group') %>% as.data.frame()
vim_tab <- bind_rows(full_AUC_res,vim_bio)

saveRDS(vim_tab,file=paste0(here::here(),"/results/SL_vim_res.RDS"))


# Fit model with just single groups of biomarkers
#loop over sets of biomarkers
SL_fits_list_individual_lab <- list()
for(i in 1:length(all_milk_components)){
  #print(i)
  SL_fits=AUC_res=NULL
  set.seed(12345)
  SL_fits = fit_SuperLearner(d,outcome="arm", covars=Xvars[(Xvars %in% all_milk_components[[i]])], slmod=sl, CV=T, family="binomial")
  AUC_res = extract_AUC(SL_fits)
  AUC_res
  
  SL_fits_list_individual_lab[[i]] <- list(SL_fits=SL_fits, AUC_res=AUC_res)
}

names(SL_fits_list_individual_lab) <- names(all_milk_components)

saveRDS(SL_fits_list_individual_lab,file=paste0(here::here(),"/data/models/SL_individual_lab_vim_fits.RDS"))

#save combined AUC tables for plotting
full_AUC_res <- full_AUC_res %>% mutate(group='all')
vim_bio_individual_lab <- rbindlist(lapply(SL_fits_list_individual_lab, function(x) x[[2]]), idcol='group') %>% as.data.frame()
vim_tab_individual_lab <- bind_rows(full_AUC_res,vim_bio_individual_lab)

saveRDS(vim_tab_individual_lab,file=paste0(here::here(),"/results/SL_individual_lab_vim_res.RDS"))
saveRDS(SL_fits_list,file=paste0(here::here(),"/results/SL_vim_fits.RDS"))

