

#https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
#https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

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
table(d$arm)
d$arm <- factor(d$arm, levels=c("Control","BEP","Nico"))
levels(d$arm)
table(d$study, d$arm)
table(is.na(d$arm))

#Check for missingness in adjustment covariates.
missing_W <- d %>% select(all_of(Wvars)) %>% summarise_all(funs(sum(is.na(.))))
missing_W      


SL.lib  = c("SL.mean","SL.glm","SL.glmnet","SL.xgboost")



res_primary <- d %>% group_by(study, visit) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=c(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit),
                     scale = FALSE,
                     bppar.type = BiocParallel::SnowParam()))
names(res_primary$res) <- paste0(res_primary$study, "-", res_primary$visit)
saveRDS(res_primary, file=paste0(here::here(),"/results/adjusted_combined_arms_primary_intervention_effects_results_unscaled.RDS"))
res_primary <- readRDS(paste0(here::here(),"/results/adjusted_combined_arms_primary_intervention_effects_results_unscaled.RDS"))

res_secondary <- d %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=c(all_milk_components$hmo, all_milk_components$protein),
                     scale = FALSE,
                     bppar.type = BiocParallel::SnowParam())))
names(res_secondary$res) <- paste0(res_secondary$study, "-", res_secondary$visit)
saveRDS(res_secondary, file=paste0(here::here(),"/results/adjusted_combined_arms_secondary_intervention_effects_results_unscaled.RDS"))
res_secondary <- readRDS(paste0(here::here(),"/results/adjusted_combined_arms_secondary_intervention_effects_results_unscaled.RDS"))

#temp simplify library
SL.lib  = c("SL.glm")

res_tertiary <- d %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T,  g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=all_milk_components$metabolomics,
                     scale = FALSE,
                     bppar.type = BiocParallel::SnowParam())))
names(res_tertiary$res) <- paste0(res_tertiary$study, "-", res_tertiary$visit)
saveRDS(res_tertiary, file=paste0(here::here(),"/results/adjusted_combined_arms_tertiary_intervention_effects_results_unscaled.RDS"))

saveRDS(list(res_primary=res_primary, 
             res_secondary=res_secondary,
             res_tertiary=res_tertiary),
        file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results_unscaled.RDS"))

