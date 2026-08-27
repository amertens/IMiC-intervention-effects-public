# =============================================================================
# src/3 visualizations/2-volcano_plots-combined-arms-traj.R
#
# Reads:  results/adjusted_intervention_effects_traj_results_clean.RDS
# Writes: figures/volcano_plot_traj_primary.png
#         figures/volcano_plot_traj_secondary.png
#         figures/volcano_plot_traj_tertiary.png
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

res <- readRDS(file=paste0(here::here(),"/results/adjusted_intervention_effects_traj_results_clean.RDS"))


res_primary <- res %>% filter(outcome_group=="primary", measure=="ATE") 
res_secondary <- res %>% filter(outcome_group=="secondary", measure=="ATE") 
res_tertiary <- res %>% filter(outcome_group=="tertiary", measure=="ATE") 



# p_primary <- plot_imic_volcano(res_primary, title="Trajectories: primary outcomes", facet_type="study")
# p_secondary <- plot_imic_volcano(res_secondary, title="Trajectories: seconday outcomes", facet_type="study")
# p_tertiary <- plot_imic_volcano(res_tertiary, title="Trajectories: tertiary outcomes", facet_type="study")
# p_primary
# p_secondary
# p_tertiary

p_primary_category <- plot_imic_volcano(res_primary, title="Trajectories: primary outcomes", facet_type="study", label_type="category")
p_secondary_category <- plot_imic_volcano(res_secondary, title="Trajectories: seconday outcomes", facet_type="study", label_type="category")
p_tertiary_category <- plot_imic_volcano(res_tertiary, title="Trajectories: tertiary outcomes", facet_type="study", label_type="category")
p_primary_category 
p_secondary_category 
p_tertiary_category 

ggsave(p_primary_category, file=paste0(here::here(),"/figures/volcano_plot_traj_primary.png"), width=18, height=10)
ggsave(p_secondary_category, file=paste0(here::here(),"/figures/volcano_plot_traj_secondary.png"), width=18, height=10)
ggsave(p_tertiary_category, file=paste0(here::here(),"/figures/volcano_plot_traj_tertiary.png"), width=18, height=10)
