

#https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
#https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))
table(d$study, d$arm)

colnames(d)

for(i in 32:ncol(d)){
    d[,i] <- d[,i]/d$fat
  
}


#Check for missingness in adjustment covariates.
missing_W <- d %>% select(all_of(Wvars)) %>% summarise_all(funs(sum(is.na(.))))
missing_W      

SL.lib  = c("SL.mean","SL.glm","SL.glmnet","SL.xgboost")
SL.lib  = c("SL.mean","SL.glm")




res_tertiary <- d %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T,  g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=all_milk_components$metabolomics,
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))

names(res_tertiary$res) <- paste0(res_tertiary$study, "-", res_tertiary$visit)

res <- extract_bioTMLE_results(res_tertiary, single_group = TRUE)


saveRDS(res, file=paste0(here::here(),"/results/fat_adjusted_metabolomics_intervention_effects_results.RDS"))

