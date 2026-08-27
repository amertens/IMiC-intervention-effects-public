
rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

#merge datasets
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))
growth <- read.csv(file=paste0(here::here(),"/data/imic_full_growth_outcomes_dataset.csv")) %>%
  filter(studyid!="CHILD") %>%
  mutate(subjid=as.character(subjid))

d <- d %>%
  distinct(studyid, subjid, subjido, arm, .keep_all = TRUE)

dim(d)
dim(growth)
df<- left_join(d, growth, by = c("subjid", "subjido", "studyid", "sex"))
dim(df)

df$sex <- factor(df$sex, levels=c("Female","Male"))


#Check for missingness in adjustment covariates.
missing_W <- df %>% select(all_of(Wvars)) %>% summarise_all(funs(sum(is.na(.))))
missing_W   

table(df$arm)

#specify analysis
# df$dummy<-1
# Wvars = c("arm","dummy")
# colnames(df)

SL.lib  = c("SL.mean","SL.glm","SL.biglasso" ,"SL.ranger","SL.xgboost")


primary_outcomes <- c("haz_birth","haz_3mo", "haz_6mo", "whz_birth","whz_3mo",  "whz_6mo")
secondary_outcomes <- c("waz_6mo", "hcaz_6mo",
                        "len_vel_z_birth_6mo", 
                        "wt_vel_z_birth_6mo", 
                        "growth_faltering_birth_6mo")
growthData <- df %>% select(all_of(primary_outcomes),all_of(secondary_outcomes))
Yvars=colnames(growthData)
Yvars

df2 <- df %>% select(all_of(Yvars), all_of(Wvars))

# #run analysis
growth_res <- df %>% group_by(study) %>%
  do(fit=run_bioTMLE(d=.,  Wvars = Wvars,
                     Yvars=Yvars,
                     scale = T))

growth_res <- df %>% group_by(study) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                         Yvars=Yvars,
                         scale = TRUE,
                         bppar.type = BiocParallel::SnowParam())))

saveRDS(growth_res, file=paste0(here::here(),"/results/growth_intervention_effects_results.RDS"))

