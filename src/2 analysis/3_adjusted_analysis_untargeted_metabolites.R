# =============================================================================
# src/2 analysis/3_adjusted_analysis_untargeted_metabolites.R
#
# Untargeted milk metabolomics: reads the per-study untargeted feature tables,
# merges them with the covariates, saves the merged dataset, and fits
# arm-stratified adjusted biotmle intervention effects (GLM-only library) for
# every feature by study and visit. clean_results.R turns the output into
# results/adjusted_intervention_effects_res_untargeted_metabolomics_clean{,_ATE}.RDS
# (arm-stratified input to Table S6). The merged dataset is also the input of
# 3_adjusted_analysis_untargeted_metabolites_combined_arms.R and of the
# untargeted PCA block in 1 data prep/4-pca-reductions.R.
#
# Inputs:  data/clean milk data/{MISAME,VITAL,ELICIT}/{M,V,E}_metabolite_s.csv
#          data/merged_analysis_datasets.RDS
# Outputs: data/clean milk data/untargeted_metabolites.RData (d_metabolomics, Yvars)
#          large-file-results/metabalomics_intervention_effects_results.RDS
# [needs restricted data]
# =============================================================================

# Method references: biotmle vignette
# https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
# and https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

library(data.table)

# fread() for the large untargeted feature tables
misame <- fread(paste0(here::here(),"/data/clean milk data/MISAME/M_metabolite_s.csv")) %>% mutate(study="Misame")
gc()
colnames(misame)[1] <- "bmid"
misame[,1:20]

vital <- fread(paste0(here::here(),"/data/clean milk data/VITAL/V_metabolite_s.csv")) %>% mutate(study="Vital")
elicit <- fread(paste0(here::here(),"/data/clean milk data/ELICIT/E_metabolite_s.csv")) %>% mutate(study="Elicit")
colnames(vital)[1] <- "bmid"
colnames(elicit)[1] <- "bmid"

gc()

vital$bmid[1:30]
elicit$bmid[1:30]

metabolomics <- bind_rows(misame, vital, elicit)

metabolomics <- metabolomics %>% select("study","bmid","visit", everything())
Yvars <- colnames(metabolomics)[-c(1:3)] 

#------------------------------------------------------------------------------
# Load imic covariates data and merge
#------------------------------------------------------------------------------

d <- readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS")) %>%
  select(study, visit, subjid, subjido, bmid,  all_of(Wvars)) 

dim(d)
dim(metabolomics)
d$visit <- as.character(d$visit)
metabolomics$visit <- as.character(metabolomics$visit)
d_metabolomics <- left_join(d, metabolomics, by=c("study","bmid","visit")) 
dim(d_metabolomics)

save(d_metabolomics, Yvars, file=paste0(here::here(),"/data/clean milk data/untargeted_metabolites.RData"))

#------------------------------------------------------------------------------
# run analysis 
#------------------------------------------------------------------------------

SL.lib  = c("SL.glm")

res <- d_metabolomics %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=Yvars,
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))
names(res$res) <- paste0(res$study, "-", res$visit)

saveRDS(res, file=paste0(here::here(),"/large-file-results/metabalomics_intervention_effects_results.RDS"))

res$res[4]
res
