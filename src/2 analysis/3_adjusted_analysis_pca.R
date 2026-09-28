# =============================================================================
# src/2 analysis/3_adjusted_analysis_pca.R
#
# Estimates adjusted intervention effects (biotmle TMLE; arms pooled into
# Control, BEP and Nico as in the combined-arms analysis) on the first principal
# component of each milk-component panel, by study and visit, using the scores
# built in 1 data prep/4-pca-reductions.R. Feeds Fig 1B.
#
# Inputs:  data/pca_analysis_datasets.RDS, metadata/milk_component.Rdata
# Outputs: results/pca_intervention_effects_results.RDS
# [needs restricted data]
# =============================================================================

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/pca_analysis_datasets.RDS"))

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

table(d$arm, is.na(d$microbiome_pca))
table(d$arm, is.na(d$untarget_metabolomics_pca))


#Check for missingness in adjustment covariates.
missing_W <- d %>% select(all_of(Wvars)) %>% summarise_all(funs(sum(is.na(.))))
missing_W      

SL.lib  = c("SL.mean","SL.glm","SL.glmnet","SL.xgboost")

set.seed(123)
res <- d %>% group_by(study, visit) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=c("macro_pca","micro_pca","bvit_pca","HMO_pca","protein_pca", "metabolomics_pca",
                             "untarget_metabolomics_pca", 
                             "microbiome_pca"),
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam()))


names(res$res) <- paste0(res$study, "-", res$visit)


res_df <- Map(extract_res, res$res) %>% rbindlist(., idcol='studytime') %>% 
  mutate(outcome_group='pca') %>% 
  as.data.frame() %>% filter(measure=="ATE")
res_df

#format 
res_df <- res_df %>% mutate(
  label_f = str_to_title(gsub("_pca","",biomarker)),
  label_f = gsub("_"," ",label_f),
  label_f = gsub("Untarget","Untargeted",label_f),
  label_f = gsub("Bvit","B-vitamins",label_f),
  label_f = gsub("Protein","Proteins",label_f),
  label_f = gsub("Micro","Micronutrients",label_f),
  label_f = gsub("Micronutrientsbiome","Microbiome",label_f),
  label_f = gsub("Macro","Macronutrients",label_f),
  label_f = gsub("Metabolomics","Targeted metabolomics",label_f),
  label_f = gsub("Hmo","HMOs",label_f),
  label_f = factor(label_f, levels=rev(c("Macronutrients", "Micronutrients", "B-vitamins",
                                     "HMOs", "Proteins",
                                     "Targeted metabolomics",  "Untargeted metabolomics", "Microbiome")))
)

# The direction of the first PC has no clear meaning without its loadings, so
# report the magnitude of the shift: flip negative effects to positive. Flipping
# negates the CI bounds and swaps them, so the lower bound stays below the upper.
res_df <- res_df %>% mutate(
  flip    = est < 0,
  cil_raw = cil,
  cil     = ifelse(flip, -ciu, cil),
  ciu     = ifelse(flip, -cil_raw, ciu),
  est     = ifelse(flip, -est, est)
) %>% select(-flip, -cil_raw)

saveRDS(res_df, file=paste0(here::here(),"/results/pca_intervention_effects_results.RDS"))
