
#------------------------------------------------------------------------------
# SuperLearner variable importance of baseline covariates compared to milk components
#------------------------------------------------------------------------------


rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
library(quadprog)
library(tidyverse)
library(ck37r)
library(data.table)

#merge datasets
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))
growth <- read.csv(file=paste0(here::here(),"/data/imic_full_growth_outcomes_dataset.csv")) %>%
  filter(studyid!="CHILD") %>%
  mutate(subjid=as.character(subjid))

d <- d %>%
  distinct(studyid, subjid, subjido, arm, .keep_all = TRUE)

dim(d)
dim(growth)
df<- left_join(d, growth, by = c("subjid", "subjido", "studyid", "sex"))
dim(df)

df$sex <- factor(df$sex, levels=c("Female","Male"))
df$contrast<-"null" #to get the function working

#Check for missingness in adjustment covariates.
missing_W <- df %>% select(all_of(Wvars)) %>% summarise_all(funs(sum(is.na(.))))
missing_W   

table(df$arm)

head(df)



stack_glm <- make_learner(
  Stack,
  lrnr_mean,
  glm_fast
)

metalearner <- make_learner(Lrnr_nnls)

sl_glm <- make_learner(Lrnr_sl,
                   learners = stack_glm,
                   metalearner = metalearner)



head(df)
df$whz_6mo
full_SL_fits = fit_SuperLearner(df, outcome="whz_6mo", covars=as.vector(unlist(Wvars[1:3])), slmod=sl_glm, CV=T, family="binomial")

res = fit_SL_fun(dat=df,  outcome = "whz_6mo", id="subjid", family="gaussian",
                       covars=as.vector(unlist(Wvars[1:3])),slmod=sl_glm,
                       CV=TRUE, folds=2)

Xvars = as.vector(unlist(all_milk_components))

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

