# =============================================================================
# src/3 visualizations/descriptive_growth_plots.R
#
# Reads:  results/imic_growth_outcomes_dataset_long.RDS
#
# Paths above were recovered from this script's syntax tree and are
# repo-relative; they resolve from the repo root via here::here().
#
# Header generated from the code itself; it makes no claim about method.
# See README.md for run order and results/ARTIFACT_MANIFEST.csv for the
# exhibit each script feeds.
# =============================================================================


rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

d <- readRDS(file=paste0(here::here(),"/results/imic_growth_outcomes_dataset_long.RDS"))
       
#line plots for haz by age grouped by studyid 
p_haz <- ggplot(d, aes(x=agedays, y=haz, group=studyid, color=studyid)) +
  geom_smooth() +
  facet_wrap(~studyid) +
  coord_cartesian(ylim=c(-3,3)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(alpha=0.25) +
  theme_minimal() + theme(legend.position = "none") +
  xlab("Age (days)") + ylab("LAZ") + ggtitle("LAZ by child age")
p_whz <- ggplot(d, aes(x=agedays, y=whz, group=studyid, color=studyid)) +
  geom_smooth() +
  facet_wrap(~studyid) +
  coord_cartesian(ylim=c(-3,3)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(alpha=0.25) +
  theme_minimal() + theme(legend.position = "none") +
  xlab("Age (days)") + ylab("WLZ") + ggtitle("WLZ by child age")

