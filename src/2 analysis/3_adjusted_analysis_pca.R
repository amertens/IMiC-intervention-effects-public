# =============================================================================
# src/2 analysis/3_adjusted_analysis_pca.R
#
# Reads:  data/pca_analysis_datasets.RDS
#         metadata/milk_component.Rdata
# Writes: results/pca_intervention_effects_results.RDS
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

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/pca_analysis_datasets.RDS"))


####XXXXXXXXXXX
# TO DO: analyze the first 3 PC's, because Liat's work shows that that explains >95% of variability
# Also predict growth (underweight, wasting, and WAZ) from PCA's and do variable importance
####XXXXXXXXXXX

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

#flip ATE if negative (as direction of first PC doesn't have much meaning without looking at weights)
res_df <- res_df %>% mutate(
  cil =ifelse(est>0, cil, -cil),
  ciu =ifelse(est>0, ciu, -ciu),
  est =ifelse(est>0, est, -est)
)

saveRDS(res_df, file=paste0(here::here(),"/results/pca_intervention_effects_results.RDS"))
