# 3b-SL_vim_individual_lab.R
# =============================================================================
# SL arm-classification VIM pipeline, split 2 of 3 (was 3-SL_vim.R).
#
# *** This is the stage that feeds the manuscript Figure 1 ML/VIM panel. ***
# (figure-data/SL_vim_plot_data.RDS, built by
#  `src/3 visualizations/5-SL_VIM_plots.R` from this script's
#  results/SL_individual_lab_vim_res.RDS, is what
#  figure-scripts/manuscript_figures/fig1-ml-vim-classifier.R actually reads.
#  The leave-one-out stage, 3c, is NOT currently on that path -- its
#  results/SL_vim_res.RDS is read by 5-SL_VIM_plots.R only in a commented-out
#  line.)
#
# For each milk-component group, fits a cross-validated SuperLearner classifier
# using ONLY that group's biomarkers as covariates (the complement of 3c's
# leave-one-OUT loop). One iteration per group in `all_milk_components`
# (currently 6: macro, micro, bvit, hmo, protein, metabolomics).
#
# Depends on the full-model fit (3a) only to append `full_AUC_res` as the
# 'all' row of the combined table -- reuses the checkpoint if present instead
# of re-fitting it.
#
# Outputs:
#   data/models/SL_individual_lab_vim_fits.RDS
#   results/SL_individual_lab_vim_res.RDS
#
# Run from the repo root:
#   Rscript "src/2 analysis/3b-SL_vim_individual_lab.R"
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

# Save the SMALL, manuscript-critical results table FIRST and unconditionally.
# The 6 SuperLearner fits above are the expensive part (hours); this table
# (`results/SL_individual_lab_vim_res.RDS`, KB-sized) is what
# `src/3 visualizations/5-SL_VIM_plots.R` -> `figure-data/SL_vim_plot_data.RDS`
# -> manuscript Fig 1 actually reads. The raw multi-GB model-object list saved
# below is NOT read by anything downstream and has previously failed to write
# (this repo lives in a OneDrive-synced folder; a ~4GB saveRDS() has died
# mid-write with "error writing to connection" and left a corrupt file on
# disk) -- losing that save must not cost us the completed fits.
full_AUC_res <- full_AUC_res %>% mutate(group='all')
vim_bio_individual_lab <- rbindlist(lapply(SL_fits_list_individual_lab, function(x) x[[2]]), idcol='group') %>% as.data.frame()
vim_tab_individual_lab <- bind_rows(full_AUC_res,vim_bio_individual_lab)

saveRDS(vim_tab_individual_lab,file=paste0(here::here(),"/results/SL_individual_lab_vim_res.RDS"))
cat("wrote results/SL_individual_lab_vim_res.RDS (manuscript-critical, saved first)\n")

# Now attempt the large raw-model-object save. Non-fatal: if it fails, the
# critical output above is already safe on disk.
big_fit_path <- paste0(here::here(),"/data/models/SL_individual_lab_vim_fits.RDS")
big_save_ok <- tryCatch({
  saveRDS(SL_fits_list_individual_lab, file=big_fit_path)
  TRUE
}, error = function(e) {
  cat("WARNING: failed to save", big_fit_path, "-", conditionMessage(e), "\n")
  cat("         (non-fatal -- results/SL_individual_lab_vim_res.RDS above is unaffected)\n")
  unlink(big_fit_path)  # remove any partial/corrupt file rather than leave it looking valid
  FALSE
})
if (big_save_ok) cat("wrote data/models/SL_individual_lab_vim_fits.RDS\n")

cat("Next: rerun `src/3 visualizations/5-SL_VIM_plots.R` to refresh figure-data/SL_vim_plot_data.RDS,\n")
cat("      then `figure-scripts/manuscript_figures/fig1-ml-vim-classifier.R` to rebuild Figure 1.\n")
