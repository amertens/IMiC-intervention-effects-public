# =============================================================================
# src/2 analysis/2_adjusted_analysis_combined_arms.R
#
# Main adjusted intervention-effect models (biotmle TMLE, adjusted for Wvars) on
# the targeted milk panels by study and visit, with trial arms pooled into
# Control, BEP and Nico: primary (macronutrients, micronutrients, B-vitamins),
# secondary (HMOs, bioactive proteins) and tertiary (targeted metabolomics).
# Vital tertiary outcomes are refit without glmnet, which fails for some of
# them. clean_results.R reads both outputs; these are the estimates behind
# Figs 2-5, Fig S3 and Tables S2-S4.
#
# Inputs:  data/merged_analysis_datasets.RDS, metadata/milk_component.Rdata
# Outputs: results/adjusted_combined_arms_{primary,secondary,tertiary}_intervention_effects_results.RDS
#          results/adjusted_combined_arms_intervention_effects_results.RDS
#          results/temp_vital_adjusted_combined_arms_tertiary_intervention_effects_results.RDS
# [needs restricted data]
# =============================================================================

# Method references: biotmle vignette
# https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
# and https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

summary(d$a.tocopherol)
summary(d$g.tocopherol)
summary(d$vitamin.A)

table(d$study, d$visit)
dim(d)
length(unique(d$subjido))

# Pool trial arms by the nutritional supplement received during lactation:
# Misame BEP/BEP and IFA/BEP -> BEP, BEP/IFA (prenatal BEP only) -> Control;
# Vital BEP arms -> BEP; Elicit Nico+Az. -> Nico, Az. (azithromycin only) -> Control.
table(d$arm)
d <- d %>% mutate(
  arm = case_when(
    arm=="Az." ~ "Control",
    arm=="BEP+ExBf+AZT" ~ "BEP",
    arm=="BEP+ExBf" ~ "BEP",
    arm=="Nico+Az." ~ "Nico",
    arm=="BEP/BEP" ~ "BEP",
    arm=="IFA/BEP" ~ "BEP",
    arm=="BEP/IFA" ~ "Control",
    arm==arm ~ arm
  )
)
d$arm <- factor(d$arm, levels=c("Control","BEP","Nico"))
levels(d$arm)
table(d$study, d$arm)
table(is.na(d$arm))

#Check for missingness in adjustment covariates.
missing_W <- d %>% select(all_of(Wvars)) %>% summarise_all(funs(sum(is.na(.))))
missing_W      

#look at number of outcomes
length(c(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit))
length(c(all_milk_components$hmo, all_milk_components$protein))
length(all_milk_components$metabolomics)

SL.lib  = c("SL.mean","SL.glm","SL.glmnet","SL.xgboost")

res_primary <- d %>% group_by(study, visit) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=c(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit),
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam()))
names(res_primary$res) <- paste0(res_primary$study, "-", res_primary$visit)
saveRDS(res_primary, file=paste0(here::here(),"/results/adjusted_combined_arms_primary_intervention_effects_results.RDS"))
res_primary <- readRDS(paste0(here::here(),"/results/adjusted_combined_arms_primary_intervention_effects_results.RDS"))

res_secondary <- d %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=c(all_milk_components$hmo, all_milk_components$protein),
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))
names(res_secondary$res) <- paste0(res_secondary$study, "-", res_secondary$visit)
saveRDS(res_secondary, file=paste0(here::here(),"/results/adjusted_combined_arms_secondary_intervention_effects_results.RDS"))
res_secondary <- readRDS(paste0(here::here(),"/results/adjusted_combined_arms_secondary_intervention_effects_results.RDS"))

res_tertiary <- d %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T,  g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=all_milk_components$metabolomics,
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))
names(res_tertiary$res) <- paste0(res_tertiary$study, "-", res_tertiary$visit)
saveRDS(res_tertiary, file=paste0(here::here(),"/results/adjusted_combined_arms_tertiary_intervention_effects_results.RDS"))
res_tertiary <- readRDS(paste0(here::here(),"/results/adjusted_combined_arms_tertiary_intervention_effects_results.RDS"))


saveRDS(list(res_primary=res_primary, 
             res_secondary=res_secondary,
             res_tertiary=res_tertiary),
        file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results.RDS"))

# glmnet fails for some Vital tertiary outcomes, so refit Vital without it.
# clean_results.R swaps these fits into the combined-arms object.
SL.lib  = c("SL.mean","SL.glm","SL.xgboost")

res_vital_tertiary <- d %>% filter(study=="Vital") %>% 
  group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T,  g_lib = SL.lib, Q_lib = SL.lib,
                         Yvars=all_milk_components$metabolomics,
                         scale = TRUE,
                         bppar.type = BiocParallel::SnowParam())))
saveRDS(res_vital_tertiary, file=paste0(here::here(),"/results/temp_vital_adjusted_combined_arms_tertiary_intervention_effects_results.RDS"))
