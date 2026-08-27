

#https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
#https://joss.theoj.org/papers/10.21105/joss.00295

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


#mutate all- transform the outcomes to a trajectory by subtracting the prior observation
unique(d$visit)
d$visit <- factor(d$visit, levels=c("1",  "2",  "3",  "5",  "40", "56"))
levels(d$visit)
d <- d %>% group_by(study, subjid, subjido) %>% arrange(visit) %>%  mutate_at(vars(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit, all_milk_components$hmo, all_milk_components$protein, all_milk_components$metabolomics), funs(. - lag(.)))

d <- d %>% filter(!(visit %in% c("1", "40")))

d %>% group_by(study, visit, arm) %>%
  summarise(mean(protein, na.rm=T))


head(df)


summary(d$protein)

#check the transformation
table(d$visit, is.na(d$protein))
d$protein[d$studyid=="MISAME-3"& d$subjid=="1001"& d$subjido=="1001"]


# SL.lib  = c("SL.mean","SL.glm","SL.biglasso")
 SL.lib  = c("SL.glm")
#SL.lib  = c("SL.mean","SL.glm","SL.glmnet","SL.xgboost")

 d %>% group_by(study, visit, arm) %>%
   summarise(mean(nr, na.rm=T))
   
 
 res_primary <- d %>% group_by(study, visit) %>% droplevels() %>%
   do(res=run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                      Yvars=c(all_milk_components$bvit),
                      scale = TRUE,
                      bppar.type = BiocParallel::SnowParam()))
 res_primary$res[1]


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

#simpler library
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

