# =============================================================================
# src/2 analysis/1-unadjusted-analysis_combined_arms.R
#
# Reads:  data/merged_analysis_datasets.RDS
#         metadata/milk_component.Rdata
# Writes: results/unadjusted_combined_arms_intervention_effects_results.RDS
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
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

# df <- d %>% filter(study=="Misame", visit ==2) %>% select(arm,iga)
# ggplot(df, aes(x=arm, y=iga)) + geom_boxplot()

d$dummy<-1
Wvars = c("arm","dummy")

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

# temp <- d %>% group_by(study, visit) %>%
#   do(res=run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=F,
#                      Yvars=c(all_milk_components$macro),
#                      scale = TRUE))


res_primary <- d %>% group_by(study, visit) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=F,
                     Yvars=c(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit),
                     scale = TRUE))
names(res_primary$res) <- paste0(res_primary$study, "-", res_primary$visit)

res_secondary <- d %>% group_by(study, visit) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars,
                     Yvars=c(all_milk_components$hmo, all_milk_components$protein),
                     scale = TRUE))
names(res_secondary$res) <- paste0(res_secondary$study, "-", res_secondary$visit)

res_tertiary <- d %>% group_by(study, visit) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars,
                     Yvars=all_milk_components$metabolomics,
                     scale = TRUE))
names(res_tertiary$res) <- paste0(res_tertiary$study, "-", res_tertiary$visit)

saveRDS(list(res_primary=res_primary, 
             res_secondary=res_secondary,
             res_tertiary=res_tertiary),
        file=paste0(here::here(),"/results/unadjusted_combined_arms_intervention_effects_results.RDS"))




