# =============================================================================
# src/2 analysis/2_adjusted_analysis_HM_trajectories.R
#
# Arm-stratified adjusted intervention effects on milk-component trajectories:
# each outcome is replaced by its change from the participant's previous visit,
# the first visit of each study (1 and 40) is dropped, and biotmle (GLM-only
# library) is fit by study and visit for the primary, secondary and tertiary
# targeted panels. No printed exhibit uses these estimates; the script is kept
# because clean_results.R reads its output and writes
# results/adjusted_intervention_effects_traj_results_clean.RDS.
#
# Inputs:  data/merged_analysis_datasets.RDS, metadata/milk_component.Rdata
# Outputs: results/adjusted_HMtraj_{primary,secondary,tertiary}_intervention_effects_results.RDS
#          results/adjusted_HMtraj_intervention_effects_results.RDS
# [needs restricted data]
# =============================================================================

# Method references: biotmle vignette
# https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
# and https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

#Check for missingness in adjustment covariates.
missing_W <- d %>% select(all_of(Wvars)) %>% summarise_all(funs(sum(is.na(.))))
missing_W      

head(d)

d %>% group_by(study, visit, arm) %>%
  summarise(mean(protein, na.rm=T))


# Transform each outcome to a trajectory: subtract the participant's previous observation.
unique(d$visit)
d$visit <- factor(d$visit, levels=c("1",  "2",  "3",  "5",  "40", "56"))
levels(d$visit)
d <- d %>% group_by(study, subjid, subjido) %>% arrange(visit) %>%  mutate_at(vars(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit, all_milk_components$hmo, all_milk_components$protein, all_milk_components$metabolomics), funs(. - lag(.)))

d <- d %>% filter(!(visit %in% c("1", "40")))

d %>% group_by(study, visit, arm) %>%
  summarise(mean(protein, na.rm=T))

summary(d$protein)

#check the transformation
table(d$visit, is.na(d$protein))

SL.lib  = c("SL.glm")

res_primary <- d %>% group_by(study, visit) %>% droplevels() %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=c(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit),
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam()))
res_primary$res[1]
names(res_primary$res) <- paste0(res_primary$study, "-", res_primary$visit)
saveRDS(res_primary, file=paste0(here::here(),"/results/adjusted_HMtraj_primary_intervention_effects_results.RDS"))

res_secondary <- d %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=c(all_milk_components$hmo, all_milk_components$protein),
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))
names(res_secondary$res) <- paste0(res_secondary$study, "-", res_secondary$visit)
saveRDS(res_secondary, file=paste0(here::here(),"/results/adjusted_HMtraj_secondary_intervention_effects_results.RDS"))

SL.lib  = c("SL.glm")

res_tertiary <- d %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T,  g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=all_milk_components$metabolomics,
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))
names(res_tertiary$res) <- paste0(res_tertiary$study, "-", res_tertiary$visit)
saveRDS(res_tertiary, file=paste0(here::here(),"/results/adjusted_HMtraj_tertiary_intervention_effects_results.RDS"))

saveRDS(list(res_primary=res_primary, 
             res_secondary=res_secondary,
             res_tertiary=res_tertiary),
        file=paste0(here::here(),"/results/adjusted_HMtraj_intervention_effects_results.RDS"))
