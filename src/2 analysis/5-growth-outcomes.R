# =============================================================================
# src/2 analysis/5-growth-outcomes.R
#
# Adjusted intervention effects (biotmle TMLE, original trial arms vs control,
# adjusted for Wvars) on child growth within the IMiC substudies, by study:
# length-for-age (HAZ) and weight-for-length (WHZ) z-scores at birth, 3 and 6
# months; weight- and head-circumference-for-age z at 6 months; length and
# weight velocity z-scores from birth to 6 months; and growth faltering.
# Feeds Fig S1.
#
# Inputs:  data/merged_analysis_datasets.RDS
#          data/imic_full_growth_outcomes_dataset.csv (from 1 data prep/2-imic_calc_growth_outcomes.R)
# Outputs: results/growth_intervention_effects_results.RDS (estimate tables only)
# [needs restricted data]
# =============================================================================

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

growth_res <- df %>% group_by(study) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                         Yvars=Yvars,
                         scale = TRUE,
                         bppar.type = BiocParallel::SnowParam())))

# Keep only the estimate table: the fitted bioTMLE object embeds participant-level data.
growth_res$res <- lapply(growth_res$res, function(r) if (inherits(r, "try-error")) r else list(res = r$res))
saveRDS(growth_res, file=paste0(here::here(),"/results/growth_intervention_effects_results.RDS"))
