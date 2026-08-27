# =============================================================================
# src/3 visualizations/2-volcano_plots.R
#
# Reads:  data/metadata/Milk_Component_Spec_IMiC_V02.csv
#         results/adjusted_intervention_effects_results.RDS
#         results/unadjusted_intervention_effects_velocity_results_blinded.RDS
# Writes: figures/figure-data/volcano_plots.RDS
#         figures/figure-data/volcano_velocity_plots.RDS
#         figures/volcano_elicit.png
#         figures/volcano_misame.png
#         figures/volcano_vital.png
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

res <- readRDS(file=paste0(here::here(),"/results/adjusted_intervention_effects_results.RDS"))
#res <- readRDS(file=paste0(here::here(),"/results/unadjusted_intervention_effects_results_blinded.RDS"))
res <- extract_bioTMLE_results(res)
table(res$sig)
table(res$sigFDR)
prop.table(table(res$sig)) * 100
prop.table(table(res$sigFDR)) * 100

labels <- read.csv(paste0(here::here(),"/data/metadata/Milk_Component_Spec_IMiC_V02.csv"))%>%
  mutate(biomarker=tolower(VarName)) %>% select(biomarker, VARLABEL, VarSubCategory) %>% 
  rename("label"="VARLABEL", "category"="VarSubCategory") %>%
  mutate(#label_f=paste0(label, " (", category, ")")
         label_f=category
         ) %>% distinct(biomarker, .keep_all = TRUE)

res <- left_join(res, labels, by="biomarker")
head(res)

res_primary <- res %>% filter(outcome_group=="primary", measure=="ATE") %>%
  group_by(studytime, biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup() %>% arrange(chi_pval_adj)
res_secondary <- res %>% filter(outcome_group=="secondary", measure=="ATE") %>%
  group_by(studytime, biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup() %>% arrange(chi_pval_adj)
res_tertiary <- res %>% filter(outcome_group=="tertiary", measure=="ATE") %>%
  group_by(studytime, biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup() %>% arrange(chi_pval_adj)


unique(res$visit)
res_last_time_point <- res %>% filter(visit %in% c("5","3","56"), measure=="ATE") %>%
  group_by(studytime, biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup() %>% arrange(chi_pval_adj)


# #Mark top 10 adjusted P-values for labeling
# res= res_tertiary
# df_top_ten <- res %>% arrange(chi_pval_adj) %>%
#   group_by(studytime, contrast) %>% 
#   summarise(chi_pval_adj = first(chi_pval_adj)) %>%
#   top_n(10, chi_pval_adj) %>%
#   filter(chi_pval_adj < 0.05) %>%
#   select(studytime, contrast)
# 
# df_marked <- df %>%
#   mutate(top_ten = ifelse(group %in% df_top_ten$group, "Top 10", "Not Top 10"))


p_primary <- plot_imic_volcano(res_primary, title="Primary outcomes", facet_type="study")
p_secondary <- plot_imic_volcano(res_secondary, title="Seconday outcomes", facet_type="study")
p_tertiary <- plot_imic_volcano(res_tertiary, title="Tertiary outcomes", facet_type="study")
p_primary
p_secondary
p_tertiary


p_final_time <- plot_imic_volcano(res_last_time_point, title="Volcano plot at final timepoint", facet_type="outcome group")
p_final_time

#-------------------------------------------------------------------------------
#  Sort plots by study
#-------------------------------------------------------------------------------

res_misame <- res %>% filter(study=="Misame", measure=="ATE") %>%
  group_by(studytime, biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup()
res_elicit <- res %>% filter(study=="Elicit", measure=="ATE") %>%
  group_by(studytime, biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup()
res_vital <- res %>% filter(study=="Vital", measure=="ATE") %>%
  group_by(studytime, biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup()


p_misame <- plot_imic_volcano(res_misame, title="Misame", facet_type="outcome group")
p_elicit <- plot_imic_volcano(res_elicit, title="Elicit", facet_type="outcome group")
p_vital <- plot_imic_volcano(res_vital, title="Vital", facet_type="outcome group")
p_misame
p_elicit
p_vital


#-------------------------------------------------------------------------------
# Facet by arm and time, split by primary, secondary, and tertiary (or combined)
#-------------------------------------------------------------------------------



res_misame_arm <- res %>% filter(study=="Misame", measure=="ATE") 
res_elicit_arm <- res %>% filter(study=="Elicit", measure=="ATE") 
res_vital_arm <- res %>% filter(study=="Vital", measure=="ATE") 


p_misame_arm_primary <- plot_imic_volcano(res_misame_arm %>% filter(outcome_group=="primary"), title="Misame", facet_type="arm")
p_elicit_arm_primary <- plot_imic_volcano(res_elicit_arm %>% filter(outcome_group=="primary"), title="Elicit", facet_type="arm")
p_vital_arm_primary <- plot_imic_volcano(res_vital_arm %>% filter(outcome_group=="primary"), title="Vital", facet_type="arm")
p_misame_arm_primary
p_elicit_arm_primary
p_vital_arm_primary

p_misame_arm_secondary <- plot_imic_volcano(res_misame_arm %>% filter(outcome_group=="secondary"), title="Misame", facet_type="arm")
p_elicit_arm_secondary <- plot_imic_volcano(res_elicit_arm %>% filter(outcome_group=="secondary"), title="Elicit", facet_type="arm")
p_vital_arm_secondary <- plot_imic_volcano(res_vital_arm %>% filter(outcome_group=="secondary"), title="Vital", facet_type="arm")
p_misame_arm_secondary
p_elicit_arm_secondary
p_vital_arm_secondary

p_misame_arm_tertiary <- plot_imic_volcano(res_misame_arm %>% filter(outcome_group=="tertiary"), title="Misame", facet_type="arm")
p_elicit_arm_tertiary <- plot_imic_volcano(res_elicit_arm %>% filter(outcome_group=="tertiary"), title="Elicit", facet_type="arm")
p_vital_arm_tertiary <- plot_imic_volcano(res_vital_arm %>% filter(outcome_group=="tertiary"), title="Vital", facet_type="arm")
p_misame_arm_tertiary
p_elicit_arm_tertiary
p_vital_arm_tertiary


#-------------------------------------------------------------------------------
# Velocity plots
#-------------------------------------------------------------------------------

res_vel <- readRDS(file=paste0(here::here(),"/results/unadjusted_intervention_effects_velocity_results_blinded.RDS"))
res_vel <- extract_bioTMLE_results(res_vel)

res_vel_primary <- res_vel %>% filter(outcome_group=="primary", measure=="ATE") %>%
  group_by( biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup() %>% arrange(chi_pval_adj)
res_vel_secondary <- res_vel %>% filter(outcome_group=="secondary", measure=="ATE") %>%
  group_by( biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup() %>% arrange(chi_pval_adj)
res_vel_tertiary <- res_vel %>% filter(outcome_group=="tertiary", measure=="ATE") %>%
  group_by( biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup() %>% arrange(chi_pval_adj)



p_vel_primary <- plot_imic_volcano(res_vel_primary, title="Primary outcomes - HM Delta", facet_type="study")
p_vel_secondary <- plot_imic_volcano(res_vel_secondary, title="Seconday outcomes - HM Delta", facet_type="study")
p_vel_tertiary <- plot_imic_volcano(res_vel_tertiary, title="Tertiary outcomes - HM Delta", facet_type="study")
p_vel_primary
p_vel_secondary
p_vel_tertiary



#-------------------------------------------------------------------------------
# Save plot objects
#-------------------------------------------------------------------------------


ggsave(p_misame, file=paste0(here::here(),"/figures/volcano_elicit.png"))
ggsave(p_elicit, file=paste0(here::here(),"/figures/volcano_misame.png"))
ggsave(p_vital, file=paste0(here::here(),"/figures/volcano_vital.png"))


saveRDS(list(p_misame, p_elicit, p_vital), file=paste0(here::here(),"/figures/figure-data/volcano_plots.RDS"))
saveRDS(list(p_vel_primary, p_vel_secondary, p_vel_tertiary), file=paste0(here::here(),"/figures/figure-data/volcano_velocity_plots.RDS"))

