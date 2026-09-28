# =============================================================================
# src/2 analysis/2_adjusted_analysis.R
#
# Arm-stratified adjusted intervention effects (biotmle TMLE, each original trial
# arm vs control, adjusted for Wvars) on the targeted milk panels, by study and
# visit: primary (macronutrients, micronutrients, B-vitamins), secondary (HMOs,
# bioactive proteins) and tertiary (targeted metabolomics). clean_results.R
# reads the output; the arm-stratified results feed Fig S4 and the stratified
# inputs of Figs 3 and 5 and Tables S3-S4.
#
# Inputs:  data/merged_analysis_datasets.RDS, metadata/milk_component.Rdata
# Outputs: results/adjusted_primary_intervention_effects_results.RDS
#          results/adjusted_intervention_effects_results.RDS
# [needs restricted data]
# =============================================================================

# Method references: biotmle vignette
# https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
# and https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))
table(d$study, d$arm)

# Inspection only: how many targeted metabolites are near-binary (<= 3 distinct
# values) and how often they are detected.
tertiary <- d %>% select(all_milk_components$metabolomics)
head(tertiary)
res_vec=rep(NA, ncol(tertiary))
for(i in 1:ncol(tertiary)){
  res_vec[i] <- length(unique(tertiary[,i]))
}

min(res_vec[res_vec>3])
286/1545 *100
tertiary_bin <- tertiary[,res_vec<=3]
prop_present <- function(x){
  sum(x>0, na.rm=T)/length(na.omit(x))
}
prop_present_vec <- sapply(tertiary_bin, prop_present)
summary(prop_present_vec)

SL.lib  = c("SL.mean","SL.glm","SL.glmnet","SL.xgboost")

res_primary <- d %>% group_by(study, visit) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=c(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit),
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam()))
names(res_primary$res) <- paste0(res_primary$study, "-", res_primary$visit)
saveRDS(res_primary, file=paste0(here::here(),"/results/adjusted_primary_intervention_effects_results.RDS"))
res_primary <- readRDS(paste0(here::here(),"/results/adjusted_primary_intervention_effects_results.RDS"))

res_secondary <- d %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=c(all_milk_components$hmo, all_milk_components$protein),
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))
names(res_secondary$res) <- paste0(res_secondary$study, "-", res_secondary$visit)


res_tertiary <- d %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T,  g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=all_milk_components$metabolomics,
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))

names(res_tertiary$res) <- paste0(res_tertiary$study, "-", res_tertiary$visit)

saveRDS(list(res_primary=res_primary, 
             res_secondary=res_secondary,
             res_tertiary=res_tertiary),
        file=paste0(here::here(),"/results/adjusted_intervention_effects_results.RDS"))
