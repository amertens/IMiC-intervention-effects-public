# =============================================================================
# src/3 visualizations/metabolics_trajectory_plots.R
#
# Reads:  data/merged_analysis_datasets.RDS
#         results/adjusted_combined_arms_intervention_effects_results_clean.RDS
#         results/adjusted_subgroup_results.RDS
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

res_full <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))
res_subgroup_full <- readRDS(file=paste0(here::here(),"/results/adjusted_subgroup_results.RDS"))

res_subgroup_full <- clean_biomarker_labels(res_subgroup_full$res_tertiary)


#find common triglyceride
res <- res_full %>% filter(measure=="ATE", study=="Vital", outcome_group=="tertiary", sigFDR==1) %>% arrange(chi_pval_adj)
res_subgroup <- res_subgroup_full %>% filter(measure=="ATE", study=="Vital", 
                                             interaction_pvalue<0.5, visit==56) 

common_biomarkers <- res$biomarker[res$biomarker %in% res_subgroup$biomarker]
common_biomarkers[1]
i=2
res %>% filter(biomarker==common_biomarkers[i])
res_subgroup %>% filter(biomarker==common_biomarkers[i])

temp =res_subgroup %>% filter(biomarker==common_biomarkers[i])

dfull<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))
head(dfull)
#d <- dfull %>% filter(study=="Vital") %>% select(study,mbmi, arm, visit, tg.16.0_36.3.)
d <- dfull %>% filter(study=="Vital") %>% select(study,mbmi, arm, visit, tg.22.6_32.0.)
colnames(d)[5] <- "biomarker"

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
d$m_underweight <- ifelse(d$mbmi <= 18.5,"underweight","normal or overweight")
d$biomarker_scaled <- scale(d$biomarker)

ggplot(d, aes(x=visit, y=biomarker_scaled, color=arm)) + 
  geom_boxplot() + theme_ki() + 
  theme(legend.position="bottom") 

ggplot(d, aes(x=visit, y=biomarker_scaled, color=arm)) + 
  geom_boxplot() +  facet_wrap(~m_underweight) + 
  theme_ki() + 
  theme(legend.position="bottom")

#Note! The BEP intervention may just be changing the compositions of triglycerides (increasing some, decreasing others in the underweight moms)

d <- dfull %>% filter(study=="Vital") %>% select(mbmi, arm, visit, starts_with("tg."))
d$m_underweight <- ifelse(d$mbmi <= 18.5,"underweight","normal or overweight")
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
#tranform to long format and summarize mean value
d_long <- d %>% gather(key=biomarker, value=value, -c(mbmi,m_underweight,arm,visit)) %>%
  group_by( biomarker) %>%
  mutate(value=scale(value)) %>%
  group_by( arm, visit, biomarker) %>%
  summarise(mean_value=mean(value, na.rm=T)) %>% 
  ungroup() %>% mutate(m_underweight="unstratified")

head(d_long)
d_long_strat <- d %>% gather(key=biomarker, value=value, -c(mbmi,m_underweight,arm,visit)) %>%
  group_by( biomarker) %>%
  mutate(value=scale(value)) %>%
  group_by( arm, m_underweight, visit, biomarker) %>%
  summarise(mean_value=mean(value, na.rm=T)) %>% 
  ungroup()

plotdf <- bind_rows(d_long, d_long_strat)

#arrange dataset by biggest differences between groups
plotdf <- plotdf %>% group_by(biomarker) %>% 
  mutate(m_underweight=factor(m_underweight, levels=c("underweight","normal or overweight","unstratified"))) %>%
  mutate(SD=sd(mean_value)) %>% 
  arrange(SD) %>% ungroup() %>%
  mutate(biomarker=factor(biomarker, levels=unique(biomarker)))

#heatmap, stratified by underweight
ggplot(plotdf, aes(x=paste0(visit, " - ", arm), y=biomarker, fill=mean_value)) + 
  geom_tile() + scale_fill_viridis_c() + 
  facet_wrap(~m_underweight, scales="free") + 
  theme_minimal() + 
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank(),
        legend.position = "none") + 
  ggtitle("Mean value of triglycerides in Vital by visit and arm")
