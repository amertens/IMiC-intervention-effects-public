# =============================================================================
# src/3 visualizations/2-volcano_plots-combined-arms.R
#
# Reads:  results/adjusted_combined_arms_and_visits_intervention_effects_results_clean.RDS
#         results/adjusted_combined_arms_intervention_effects_results_clean.RDS
# Writes: figures/figure-data/volcano_plots_pooled_arms.RDS
#         figures/volcano_plot_key_timepoint.jpeg
#         figures/volcano_plot_tertiary.png
#         figures/volcano_plot_tertiary_cat_lab.png
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

res_pooled <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_and_visits_intervention_effects_results_clean.RDS"))
res <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))

# The "and_visits" pooled file ships with chi_pval / chi_pval_adj only.
# `plot_imic_volcano()` looks for pval / pval_adj — alias them so the pooled
# plots (lines 82-84) work with the same function signature.
if (!"pval"     %in% colnames(res_pooled)) res_pooled$pval     <- res_pooled$chi_pval
if (!"pval_adj" %in% colnames(res_pooled)) res_pooled$pval_adj <- res_pooled$chi_pval_adj


unique(res$label_f[grepl("holine",res$label_f)])


res_primary <- res %>% filter(outcome_group=="primary", measure=="ATE") 
res_secondary <- res %>% filter(outcome_group=="secondary", measure=="ATE") 
res_tertiary <- res %>% filter(outcome_group=="tertiary", measure=="ATE") 

res_last_time_point <- res %>% filter(visit %in% c("2 mo.","3-4 mo.","5 mo."), measure=="ATE") %>%
  group_by(studytime, biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup() %>% arrange(chi_pval_adj)

res_key_time_point <- res %>% filter(visit %in% c("1.5 mo.","3-4 mo.","5 mo."), measure=="ATE") %>%
  group_by(studytime, biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup() %>% arrange(chi_pval_adj)


p_primary <- plot_imic_volcano(res_primary, title="Primary outcomes", facet_type="study")
p_secondary <- plot_imic_volcano(res_secondary, title="Seconday outcomes", facet_type="study")
p_tertiary <- plot_imic_volcano(res_tertiary, title="Tertiary outcomes", facet_type="study")
p_primary
p_secondary

p_primary_category <- plot_imic_volcano(res_primary, title="Primary outcomes", facet_type="study", label_type="category")
p_secondary_category <- plot_imic_volcano(res_secondary, title="Seconday outcomes", facet_type="study", label_type="category")
p_tertiary_category <- plot_imic_volcano(res_tertiary, title="Tertiary outcomes", facet_type="study", label_type="category")
p_primary_category 
p_secondary_category 
p_tertiary_category 

overlap_n=12
p_tertiary_plot1 <- plot_imic_volcano(res_tertiary %>% filter(study=="Elicit"), overlap_n=overlap_n)
p_tertiary_plot2 <- plot_imic_volcano(res_tertiary %>% filter(study=="Vital"), overlap_n=overlap_n)
p_tertiary_plot3 <- plot_imic_volcano(res_tertiary %>% filter(study=="Misame"), overlap_n=overlap_n)
p_tertiary_plot <- cowplot::plot_grid(p_tertiary_plot3,p_tertiary_plot2,p_tertiary_plot1, ncol=1, rel_heights=c(1.5,1.2,1))
ggsave(p_tertiary_plot, file=paste0(here::here(),"/figures/volcano_plot_tertiary.png"), width=18, height=10)
ggsave(p_tertiary_plot, file=paste0(here::here(),"/figures/volcano_plot_tertiary.png"), width=10, height=6)

overlap_n=20
p_tertiary_plot1 <- plot_imic_volcano(res_tertiary %>% filter(study=="Elicit"),label_type="category", overlap_n=overlap_n)
p_tertiary_plot2 <- plot_imic_volcano(res_tertiary %>% filter(study=="Vital"),label_type="category", overlap_n=overlap_n)
p_tertiary_plot3 <- plot_imic_volcano(res_tertiary %>% filter(study=="Misame"),label_type="category", overlap_n=overlap_n)
p_tertiary_plot <- cowplot::plot_grid(p_tertiary_plot3,p_tertiary_plot2,p_tertiary_plot1, ncol=1, rel_heights=c(1.5,1.2,1))
ggsave(p_tertiary_plot, file=paste0(here::here(),"/figures/volcano_plot_tertiary_cat_lab.png"), width=18, height=10)


#other plot organizations
p_final_time <- plot_imic_volcano(res_last_time_point, title="Volcano plot at final timepoint", facet_type="outcome group")
p_final_time

res_key_time_point$label_f <- res_key_time_point$category
res_key_time_point$label_f <- gsub("Triglycerides","Tg.",res_key_time_point$label_f)
res_key_time_point$label_f <- gsub("Diglycerides","Dg.",res_key_time_point$label_f)
res_key_time_point$label_f <- gsub("Amino Acid Related","Amino Acid",res_key_time_point$label_f)
res_key_time_point$label_f <- gsub("B3 or related","B3",res_key_time_point$label_f)
p_key_time <- plot_imic_volcano(res_key_time_point, title="", facet_type="outcome group")
p_key_time
ggsave(p_key_time, file=paste0(here::here(),"/figures/volcano_plot_key_timepoint.jpeg"), width=9, height=7)



res_pooled_primary <- res_pooled %>% filter(outcome_group=="primary", measure=="ATE") 
res_pooled_secondary <- res_pooled %>% filter(outcome_group=="secondary", measure=="ATE") 
res_pooled_tertiary <- res_pooled %>% filter(outcome_group=="tertiary", measure=="ATE") 


# res_pooled_primary <- res_pooled %>% filter(outcome_group=="primary", measure=="ATE") 
# res_pooled_secondary <- res_pooled %>% filter(outcome_group=="secondary", measure=="ATE") 
# res_pooled_tertiary <- res_pooled %>% filter(outcome_group=="tertiary", measure=="ATE") 
# 
# p_primary_pooled <- plot_imic_volcano(res_pooled_primary, title="Primary outcomes", facet_type="study")
# p_secondary_pooled <- plot_imic_volcano(res_pooled_secondary, title="Seconday outcomes", facet_type="study")
# p_tertiary_pooled <- plot_imic_volcano(res_pooled_tertiary, title="Tertiary outcomes", facet_type="study")
# p_primary_pooled
# p_secondary_pooled
# p_tertiary_pooled
p_primary_pooled <- plot_imic_volcano(res_pooled_primary, title="Primary outcomes", facet_type="study")
p_secondary_pooled <- plot_imic_volcano(res_pooled_secondary, title="Seconday outcomes", facet_type="study")
p_tertiary_pooled <- plot_imic_volcano(res_pooled_tertiary, title="Tertiary outcomes", facet_type="study")
p_primary_pooled
p_secondary_pooled
p_tertiary_pooled

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
# Save plot objects
#-------------------------------------------------------------------------------


saveRDS(list(p_misame, p_elicit, p_vital), file=paste0(here::here(),"/figures/figure-data/volcano_plots_pooled_arms.RDS"))

