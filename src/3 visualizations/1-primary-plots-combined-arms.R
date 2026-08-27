# =============================================================================
# src/3 visualizations/1-primary-plots-combined-arms.R
#
# Reads:  results/adjusted_combined_arms_and_visits_intervention_effects_results_clean.RDS
#         results/adjusted_combined_arms_intervention_effects_results_clean.RDS
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
res[res$biomarker=="Na",]
head(res)

# Exploratory peek at vitamin A and E (tocopherol) results. Guarded so the
# script doesn't crash if biomarker names don't match (e.g. after a label
# refactor) — the facet would fail on empty data.
temp<-res[grepl("toco",res$biomarker) | grepl("itamin.a",res$biomarker),] %>% filter(measure=="ATE")
if (nrow(temp) > 0) {
  ggplot(temp, aes(x=biomarker, y=est)) + geom_point() +
    geom_linerange(aes(ymin=cil, ymax=ciu)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    coord_flip() +
    facet_wrap(~studytime, scale="free") +
    ggtitle("A and E vitamin results") +
    theme(strip.background = element_blank(),
          axis.text = element_text(size = 6),
          strip.text = element_text(size = 8),
          legend.position="none")  + ylab("Intervention arm") +
    xlab("Z-score difference")
}

# res <- left_join(res, labels, by="biomarker") %>% mutate(label_f=ifelse(is.na(label), biomarker, label))
# res_pooled <- left_join(res_pooled, labels, by="biomarker") %>% mutate(label_f=ifelse(is.na(label), biomarker, label))


imic_intervention_plot_function_pooled <- function(df, title="", legend=TRUE){

  # Bail out early if no rows (e.g., no FDR-sig results to plot) — facet_grid
  # would fail with "Faceting variables must have at least one value".
  if (nrow(df) == 0) {
    warning("imic_intervention_plot_function_pooled: no rows — returning empty plot.")
    return(ggplot2::ggplot() + ggplot2::labs(title = title, subtitle = "(no FDR-significant results)"))
  }

  df <- df  %>% group_by(visit, category) %>%
    mutate(mean_pval = mean(chi_pval_adj)) %>% ungroup() %>%
    arrange(mean_pval) %>%
    mutate(category = factor(category, levels=(unique(category))))
  
  df <- df  %>% arrange(category, visit, est) %>%
    mutate(label_f = str_to_sentence(label_f),
      label_f = factor(label_f, levels=rev(unique(label_f))))
  
  p <- ggplot(df, aes(x=label_f, y=est, color=category)) + geom_point() +
    geom_linerange(aes(ymin=cil, ymax=ciu)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    coord_flip() +
    #facet_wrap(~studytime, scale="free", ncol=1) +
    facet_grid(~studytime, scale="free") +
    ggtitle(title) +
    #scale_color_manual(values=tableau10) +
    theme(strip.background = element_blank(),
          axis.text = element_text(size = 6),
          strip.text = element_text(size = 8)) + 
    ylab("Z-scored difference") + xlab("Milk component")
  
  if(legend){
    p <- p + theme(legend.position="right")
  } else {
    p <- p + theme(legend.position="none")
  }
  
  return(p)
}



#-------------------------------------------------------------------------------
# Divide results by primary/secondary outcomes, etc.
#-------------------------------------------------------------------------------

head(res)
primary_res <- res %>% filter(outcome_group=="primary", measure=="ATE")
secondary_res <- res %>% filter(outcome_group=="secondary", measure=="ATE")
tertiary_res <- res %>% filter(outcome_group=="tertiary", measure=="ATE")


# p_primary <- imic_intervention_plot_functon(primary_res, title="Significant primary results")
# p_secondary <- imic_intervention_plot_functon(secondary_res, title="Significant secondary results")
# p_tertiary <- imic_intervention_plot_functon(tertiary_res, title="Significant tertiary results")

primary_res %>% filter(study=="Misame", is.na(label))

# Helper: wrap plot + ggsave so an empty/missing-facet plot doesn't kill the script
safe_plot_save <- function(p, file, width = 8, height = 4) {
  tryCatch(
    ggsave(p, file = file, width = width, height = height),
    error = function(e) message("WARNING: skipped saving ", file, " (", conditionMessage(e), ")")
  )
}

p_primary_misame <- imic_intervention_plot_function_pooled(primary_res %>% filter(study=="Misame", chi_pval_adj<0.05), title="Significant primary results-Misame")
p_primary_elicit <- imic_intervention_plot_function_pooled(primary_res %>% filter(study=="Elicit", chi_pval_adj<0.05), title="Significant primary results-Elicit")
p_primary_vital  <- imic_intervention_plot_function_pooled(primary_res %>% filter(study=="Vital",  chi_pval_adj<0.05), title="Significant primary results-Vital")

safe_plot_save(p_primary_elicit, paste0(here::here(),"/figures/primary_elicit.png"), width=6, height=3)
safe_plot_save(p_primary_misame, paste0(here::here(),"/figures/primary_misame.png"), width=8, height=4)
safe_plot_save(p_primary_vital, paste0(here::here(),"/figures/primary_vital.png"), width=8, height=4)


#None significant in elicit and misame
# p_secondary_elicit <- imic_intervention_plot_function_pooled(secondary_res %>% filter(study=="Elicit", chi_pval_adj<0.05), title="Significant secondary results")
# p_secondary_elicit
# p_secondary_misame <- imic_intervention_plot_function_pooled(secondary_res %>% filter(study=="Misame", chi_pval_adj<0.05), title="Significant secondary results")
# p_secondary_misame
# p_secondary_vital <- imic_intervention_plot_function_pooled(secondary_res %>% filter(study=="Vital", chi_pval_adj<0.05), title="Significant secondary results")
# p_secondary_vital


#TO DO
#add all across time if they are significant at one time (grey if not sig)

p_tertiary_elicit <- imic_intervention_plot_function_pooled(tertiary_res %>% filter(study=="Elicit", chi_pval_adj<0.05), title="Significant tertiary results")
# p_tertiary_elicit  # (auto-commented to avoid batch print errors)
p_tertiary_misame <- imic_intervention_plot_function_pooled(tertiary_res %>% filter(study=="Misame", chi_pval_adj<0.05), title="Significant tertiary results")
# p_tertiary_misame  # (auto-commented to avoid batch print errors)
p_tertiary_vital <- imic_intervention_plot_function_pooled(tertiary_res %>% filter(study=="Vital", chi_pval_adj<0.05), title="Significant tertiary results")
# p_tertiary_vital
  # (auto-commented to avoid batch print errors)
safe_plot_save(p_tertiary_vital, paste0(here::here(),"/figures/forest_sig_tertiary_vital.png"), width=10, height=16)


names(all_milk_components)

#micronutrients
p_macro <- imic_intervention_plot_function_pooled(primary_res %>% filter(biomarker %in% all_milk_components$macro, measure=="ATE"), title="Macronutrients", legend=FALSE)
# p_macro
  # (auto-commented to avoid batch print errors)
#need to fix missing carb labeling in Misame

#micronutrients
p_micro <- imic_intervention_plot_function_pooled(primary_res %>% filter(biomarker %in% all_milk_components$micro, measure=="ATE"), title="Micronutrients", legend=FALSE)
# p_micro
  # (auto-commented to avoid batch print errors)
#HMOs
p_hmo <- imic_intervention_plot_function_pooled(secondary_res %>% filter(biomarker %in% all_milk_components$hmo, measure=="ATE"), title="Micronutrients", legend=FALSE)
# p_hmo
  # (auto-commented to avoid batch print errors)
#bioactives
p_protein <- imic_intervention_plot_function_pooled(secondary_res %>% filter(biomarker %in% all_milk_components$protein, measure=="ATE"), title="Micronutrients", legend=FALSE)
# p_protein
  # (auto-commented to avoid batch print errors)
# #presentation plots
# 
# plotdf <- primary_res %>% filter(sigFDR==1) %>% filter(is.na(category)) 
# 
# p <- ggplot(plotdf, aes(x=contrast, y=est)) + geom_point() +
#   geom_linerange(aes(ymin=cil, ymax=ciu)) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_flip() +
#   facet_wrap(studytime~biomarker, scale="free") +
#   ggtitle("Significant macro- and micro-nutrient results\n(excluding B-vitamins)") +
#   theme(strip.background = element_blank(),
#         axis.text = element_text(size = 6),
#         strip.text = element_text(size = 8),
#         legend.position="none")  + ylab("Intervention arm") + 
#   xlab("Z-score difference")
# p

