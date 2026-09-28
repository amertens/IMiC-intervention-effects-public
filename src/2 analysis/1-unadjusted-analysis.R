# =============================================================================
# src/2 analysis/1-unadjusted-analysis.R
#
# Unadjusted, arm-stratified intervention effects (biotmle with only the arm and
# a constant dummy as covariates) on the primary, secondary and tertiary
# targeted milk panels, by study and visit. No printed exhibit uses these
# estimates directly; the script is kept because clean_results.R reads its
# output and writes results/unadjusted_intervention_effects_results_clean.RDS.
#
# Inputs:  data/merged_analysis_datasets.RDS, metadata/milk_component.Rdata
# Outputs: results/unadjusted_intervention_effects_results.RDS
# [needs restricted data]
# =============================================================================

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

# Unadjusted: the covariate set is the arm plus a constant.
d$dummy<-1
Wvars = c("arm","dummy")

res_primary <- d %>% group_by(study, visit) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=F,
                     Yvars=c(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit),
                     scale = TRUE))
names(res_primary$res) <- paste0(res_primary$study, "-", res_primary$visit)

res_secondary <- d %>% group_by(study, visit) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars,
                     Yvars=c(all_milk_components$hmo, all_milk_components$protein),
                     scale = TRUE))
names(res_secondary$res) <- paste0(res_secondary$study, "-", res_secondary$visit)

res_tertiary <- d %>% group_by(study, visit) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars,
                     Yvars=all_milk_components$metabolomics,
                     scale = TRUE))
names(res_tertiary$res) <- paste0(res_tertiary$study, "-", res_tertiary$visit)

saveRDS(list(res_primary=res_primary, 
             res_secondary=res_secondary,
             res_tertiary=res_tertiary),
        file=paste0(here::here(),"/results/unadjusted_intervention_effects_results.RDS"))
