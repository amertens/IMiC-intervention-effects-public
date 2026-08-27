# 3c-SL_vim_leave_one_out.R
# =============================================================================
# SL arm-classification VIM pipeline, split 3 of 3 (was 3-SL_vim.R).
#
# NOT currently on the path to any manuscript figure -- `5-SL_VIM_plots.R`
# reads results/SL_individual_lab_vim_res.RDS (3b's output); the read of this
# script's results/SL_vim_res.RDS is commented out there. Kept/split out for
# completeness and in case that's revisited.
#
# For each milk-component group, fits a cross-validated SuperLearner classifier
# using every OTHER group's biomarkers as covariates (leave-that-group-out;
# the complement of 3b's single-group loop). One iteration per group in
# `all_milk_components` (currently 6: macro, micro, bvit, hmo, protein,
# metabolomics).
#
# Depends on the full-model fit (3a) only to append `full_AUC_res` as the
# 'all' row of the combined table -- reuses the checkpoint if present instead
# of re-fitting it.
#
# Outputs:
#   data/models/SL_vim_fits.RDS
#   results/SL_vim_res.RDS
#   results/SL_vim_fits.RDS   (duplicate save, preserved from the original script)
#
# Run from the repo root:
#   Rscript "src/2 analysis/3c-SL_vim_leave_one_out.R"
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

#### full-model AUC row (reuse 3a's checkpoint if present) --------------------
full_fit_path <- paste0(here::here(),"/data/models/SL_vim_full_fit.RDS")
if (file.exists(full_fit_path)) {
  cat("Reusing existing full-model checkpoint:", full_fit_path, "\n")
  full_AUC_res <- readRDS(full_fit_path)$full_AUC_res
} else {
  cat("No full-model checkpoint found -- fitting it now (see 3a-SL_vim_full_model.R).\n")
  set.seed(12345)
  full_SL_fits = fit_SuperLearner(d,outcome="arm", covars=as.vector(unlist(all_milk_components)), slmod=sl, CV=T, family="binomial")
  full_AUC_res = extract_AUC(full_SL_fits)
  saveRDS(list(full_AUC_res=full_AUC_res, full_SL_fits=full_SL_fits), file=full_fit_path)
}

#### leave-one-group-out fits --------------------------------------------------
SL_fits_list <- list()
for(i in 1:length(all_milk_components)){
  cat(sprintf("[%d/%d] leaving out group: %s\n", i, length(all_milk_components), names(all_milk_components)[i]))
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
saveRDS(SL_fits_list,file=paste0(here::here(),"/results/SL_vim_fits.RDS"))

cat("wrote data/models/SL_vim_fits.RDS\n")
cat("wrote results/SL_vim_res.RDS\n")
cat("wrote results/SL_vim_fits.RDS\n")
