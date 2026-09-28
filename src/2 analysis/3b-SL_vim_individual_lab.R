# =============================================================================
# src/2 analysis/3b-SL_vim_individual_lab.R
#
# Arm-classification variable importance for Fig 1A. For each trial contrast
# (one intervention arm vs its control, by study and visit), fits a
# cross-validated sl3 SuperLearner classifier of arm from the milk components
# and records the cross-validated AUC: once with all components (the "all"
# row) and once per component group in all_milk_components (macronutrients,
# micronutrients, B-vitamins, HMOs, bioactive proteins, targeted
# metabolomics). The all-components fit is read from data/models/SL_vim_full_fit.RDS
# when that file exists; otherwise it is fit here and saved there.
# src/3 visualizations/5-SL_VIM_plots.R turns the AUC table into
# figure-data/SL_vim_plot_data.RDS, which fig1-ml-vim-classifier.R plots.
#
# Inputs:  data/merged_analysis_datasets.RDS, metadata/milk_component.Rdata
#          data/models/SL_vim_full_fit.RDS (optional; reused if present)
# Outputs: results/SL_individual_lab_vim_res.RDS
#          data/models/SL_vim_full_fit.RDS (when fit here)
#          data/models/SL_individual_lab_vim_fits.RDS (fitted models; not read downstream)
# Run from the repo root: Rscript "src/2 analysis/3b-SL_vim_individual_lab.R"
# [needs restricted data]
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

# Stack one copy of the data per contrast (intervention arm + its control) so
# fit_SuperLearner() can group by contrast.
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

#### all-components AUC row (reuse the saved fit if present) ------------------
full_fit_path <- paste0(here::here(),"/data/models/SL_vim_full_fit.RDS")
if (file.exists(full_fit_path)) {
  cat("Reusing existing full-model checkpoint:", full_fit_path, "\n")
  full_AUC_res <- readRDS(full_fit_path)$full_AUC_res
} else {
  cat("No full-model checkpoint found -- fitting it now.\n")
  set.seed(12345)
  full_SL_fits = fit_SuperLearner(d,outcome="arm", covars=as.vector(unlist(all_milk_components)), slmod=sl, CV=T, family="binomial")
  full_AUC_res = extract_AUC(full_SL_fits)
  saveRDS(list(full_AUC_res=full_AUC_res, full_SL_fits=full_SL_fits), file=full_fit_path)
}

#### single-biomarker-group fits (individual "lab") ---------------------------
SL_fits_list_individual_lab <- list()
for(i in 1:length(all_milk_components)){
  cat(sprintf("[%d/%d] fitting group: %s\n", i, length(all_milk_components), names(all_milk_components)[i]))
  SL_fits=AUC_res=NULL
  set.seed(12345)
  SL_fits = fit_SuperLearner(d,outcome="arm", covars=Xvars[(Xvars %in% all_milk_components[[i]])], slmod=sl, CV=T, family="binomial")
  AUC_res = extract_AUC(SL_fits)
  AUC_res

  SL_fits_list_individual_lab[[i]] <- list(SL_fits=SL_fits, AUC_res=AUC_res)
}

names(SL_fits_list_individual_lab) <- names(all_milk_components)

# Save the small AUC table (the Fig 1A input) before the multi-GB model list.
# The model list is not read downstream, and a failed write of that large file
# should not lose the completed fits.
full_AUC_res <- full_AUC_res %>% mutate(group='all')
vim_bio_individual_lab <- rbindlist(lapply(SL_fits_list_individual_lab, function(x) x[[2]]), idcol='group') %>% as.data.frame()
vim_tab_individual_lab <- bind_rows(full_AUC_res,vim_bio_individual_lab)

saveRDS(vim_tab_individual_lab,file=paste0(here::here(),"/results/SL_individual_lab_vim_res.RDS"))
cat("wrote results/SL_individual_lab_vim_res.RDS\n")

# Save the fitted models. Non-fatal: a failure here leaves the AUC table intact.
big_fit_path <- paste0(here::here(),"/data/models/SL_individual_lab_vim_fits.RDS")
big_save_ok <- tryCatch({
  saveRDS(SL_fits_list_individual_lab, file=big_fit_path)
  TRUE
}, error = function(e) {
  cat("Warning: failed to save", big_fit_path, "-", conditionMessage(e), "\n")
  cat("         (non-fatal; results/SL_individual_lab_vim_res.RDS is unaffected)\n")
  unlink(big_fit_path)  # remove any partial/corrupt file rather than leave it looking valid
  FALSE
})
if (big_save_ok) cat("wrote data/models/SL_individual_lab_vim_fits.RDS\n")

cat("Next: rerun `src/3 visualizations/5-SL_VIM_plots.R` to refresh figure-data/SL_vim_plot_data.RDS,\n")
cat("      then `figure-scripts/manuscript_figures/fig1-ml-vim-classifier.R` to rebuild Fig 1.\n")
