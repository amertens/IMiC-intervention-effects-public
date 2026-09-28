# =============================================================================
# src/2 analysis/clean_results.R
#
# Collects the raw bioTMLE result objects written by the milk analysis scripts
# (1-unadjusted, 2_adjusted*, 3_adjusted*), flattens each into a tidy data frame
# with extract_bioTMLE_results(), and saves the clean tables the figures and
# downstream scripts read. Its outputs feed Fig 2 (subsetted results/primary_*),
# Fig S3 (subsetted results/secondary_*), Figs 3 and 5 and Tables S2-S4
# (combined_intervention_effects_results_*_arms.RDS), Fig 4 and the blood
# scripts (adjusted_combined_arms_..._results_clean.RDS), Fig S4
# (adjusted_intervention_effects_results_clean.RDS), Table S1 (unscaled clean
# RDS), Fig 6 and Tables S5-S8 (untargeted and proteomics clean RDS). It reads
# the output of every 1-/2_/3_ milk analysis script, including the unadjusted
# and trajectory variants that no printed exhibit uses, so all of those scripts
# must run first; run the blood scripts (12 onward) after it.
#
# Inputs:  results/adjusted_combined_arms_intervention_effects_results.RDS
#          results/temp_vital_adjusted_combined_arms_tertiary_intervention_effects_results.RDS
#          results/adjusted_intervention_effects_results.RDS
#          results/unadjusted_intervention_effects_results.RDS
#          results/adjusted_{combined_arms_,}intervention_effects_results_unscaled.RDS
#          results/adjusted_{combined_arms_,}HMtraj_intervention_effects_results.RDS
#          results/microbiome_intervention_effects_results{,_arm_strat}.RDS
#          results/proteomics_intervention_effects_results{,_combined_arms}.RDS
#          large-file-results/metabalomics_intervention_effects_results{,_combined_arms}.RDS
#          metadata/milk_component.Rdata
# Outputs: results/{adjusted_combined_arms_,adjusted_,unadjusted_}intervention_effects_results_clean.RDS
#          results/adjusted_combined_arms_intervention_effects_results_clean.csv
#          results/adjusted_{combined_arms_,}intervention_effects_{traj,unscaled}_results_clean.RDS
#          results/adjusted_combined_arms_intervention_effects_{proteomics,untargeted}_results_clean{,_ATE,_MN}.RDS
#          results/adjusted_intervention_effects_res_untargeted_metabolomics_clean{,_ATE,_MN}.RDS
#          results/combined_intervention_effects_results_{combined,stratified}_arms.RDS
#          results/subsetted results/{primary_macro,primary_micro,primary_bvit,secondary_hmo,
#            secondary_bioactives,tertiary_targeted_metabolomics}{,_arm_strat}.csv
# [needs restricted data]
# =============================================================================

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))

res <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results.RDS"))
# Vital tertiary (targeted metabolomics) results come from the rerun without
# glmnet in 2_adjusted_analysis_combined_arms.R, because glmnet fails for some
# Vital outcomes; slot them into the combined-arms object.
res_vital_tertiary <- readRDS(file=paste0(here::here(),"/results/temp_vital_adjusted_combined_arms_tertiary_intervention_effects_results.RDS"))
res$res_tertiary$res$`Vital-40`$res <- res_vital_tertiary$res[[1]]$res
res$res_tertiary$res$`Vital-56`$res <- res_vital_tertiary$res[[2]]$res

res_arm_strat <- readRDS(file=paste0(here::here(),"/results/adjusted_intervention_effects_results.RDS"))
res_unadjusted_arm_strat <- readRDS(file=paste0(here::here(),"/results/unadjusted_intervention_effects_results.RDS"))
res_unscaled_arm_strat <- readRDS(file=paste0(here::here(),"/results/adjusted_intervention_effects_results_unscaled.RDS"))
res_unscaled <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results_unscaled.RDS"))
res_traj_arm_strat <- readRDS(file=paste0(here::here(),"/results/adjusted_HMtraj_intervention_effects_results.RDS"))
res_traj <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_HMtraj_intervention_effects_results.RDS"))

res_microbiome <- readRDS(paste0(here::here(),"/results/microbiome_intervention_effects_results.RDS"))
res_microbiome_arm_strat <- readRDS(paste0(here::here(),"/results/microbiome_intervention_effects_results_arm_strat.RDS"))
res_protien <- readRDS(paste0(here::here(),"/results/proteomics_intervention_effects_results_combined_arms.RDS"))
res_protien_arm_strat <- readRDS(paste0(here::here(),"/results/proteomics_intervention_effects_results.RDS"))

res_untargeted_metabolomics <- readRDS(paste0(here::here(),"/large-file-results/metabalomics_intervention_effects_results_combined_arms.RDS"))
res_untargeted_metabolomics_arm_strat <- readRDS(paste0(here::here(),"/large-file-results/metabalomics_intervention_effects_results.RDS"))

res <- extract_bioTMLE_results(res)
res_arm_strat <- extract_bioTMLE_results(res_arm_strat)
res_unadjusted_arm_strat <- extract_bioTMLE_results(res_unadjusted_arm_strat)
res_unscaled <- extract_bioTMLE_results(res_unscaled)
res_unscaled_arm_strat <- extract_bioTMLE_results(res_unscaled_arm_strat)
res_traj <- extract_bioTMLE_results(res_traj)
res_traj_arm_strat <- extract_bioTMLE_results(res_traj_arm_strat)
res_untargeted_metabolomics <- extract_bioTMLE_results(res_untargeted_metabolomics, single_group = TRUE)
res_untargeted_metabolomics_arm_strat <- extract_bioTMLE_results(res_untargeted_metabolomics_arm_strat, single_group = TRUE)

res_microbiome <- extract_bioTMLE_results(res_microbiome, single_group = TRUE)
res_microbiome_arm_strat <- extract_bioTMLE_results(res_microbiome_arm_strat, single_group = TRUE)

res_protien <- extract_bioTMLE_results(res_protien, single_group = TRUE)
res_protien_arm_strat <- extract_bioTMLE_results(res_protien_arm_strat, single_group = TRUE)
res_protien$label <- res_protien_arm_strat$label <- res_protien$description <- res_protien_arm_strat$description <- res_protien$category  <- res_protien_arm_strat$category  <- res_protien$category_raw <- res_protien_arm_strat$category_raw <- "Proteomics"

saveRDS(res,file=paste0(here::here(), "/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))
saveRDS(res_arm_strat,file=paste0(here::here(), "/results/adjusted_intervention_effects_results_clean.RDS"))
saveRDS(res_unadjusted_arm_strat,file=paste0(here::here(), "/results/unadjusted_intervention_effects_results_clean.RDS"))
saveRDS(res_traj,file=paste0(here::here(), "/results/adjusted_combined_arms_intervention_effects_traj_results_clean.RDS"))
saveRDS(res_traj_arm_strat,file=paste0(here::here(), "/results/adjusted_intervention_effects_traj_results_clean.RDS"))
saveRDS(res_unscaled,file=paste0(here::here(), "/results/adjusted_combined_arms_intervention_effects_unscaled_results_clean.RDS"))
saveRDS(res_unscaled_arm_strat,file=paste0(here::here(), "/results/adjusted_intervention_effects_unscaled_results_clean.RDS"))

saveRDS(res_protien,file=paste0(here::here(), "/results/adjusted_combined_arms_intervention_effects_proteomics_results_clean.RDS"))
saveRDS(res_untargeted_metabolomics,file=paste0(here::here(), "/results/adjusted_combined_arms_intervention_effects_untargeted_results_clean.RDS"))
saveRDS(res_untargeted_metabolomics_arm_strat,file=paste0(here::here(), "/results/adjusted_intervention_effects_res_untargeted_metabolomics_clean.RDS"))

# The high-dimensional clean frames mix ATE rows (intervention effects) with MN
# rows (arm means). Counts of significant intervention effects must use the ATE
# rows only, so also save *_ATE.RDS and *_MN.RDS companions that downstream
# scripts read directly.
split_ate_mn <- function(obj, stem) {
  if (!"measure" %in% names(obj)) { warning("split_ate_mn: no 'measure' column in ", stem); return(invisible()) }
  ate <- obj[obj$measure == "ATE", , drop = FALSE]
  mn  <- obj[obj$measure == "MN",  , drop = FALSE]
  saveRDS(ate, paste0(here::here(), "/results/", stem, "_ATE.RDS"))
  saveRDS(mn,  paste0(here::here(), "/results/", stem, "_MN.RDS"))
  cat(sprintf("[split ATE/MN] %s: %d ATE + %d MN rows\n", stem, nrow(ate), nrow(mn)))
}
split_ate_mn(res_protien,                           "adjusted_combined_arms_intervention_effects_proteomics_results_clean")
split_ate_mn(res_untargeted_metabolomics,           "adjusted_combined_arms_intervention_effects_untargeted_results_clean")
split_ate_mn(res_untargeted_metabolomics_arm_strat, "adjusted_intervention_effects_res_untargeted_metabolomics_clean")

write.csv(res,file=paste0(here::here(),
                          "/results/adjusted_combined_arms_intervention_effects_results_clean.csv"))

# Combined result tables across all milk modalities (targeted, untargeted
# metabolomics, microbiome, proteomics), arm-stratified and combined-arms.
res_combined_arm_strat <- bind_rows(res_arm_strat, res_untargeted_metabolomics_arm_strat, res_microbiome_arm_strat, res_protien_arm_strat)
saveRDS(res_combined_arm_strat,file=paste0(here::here(),"/results/combined_intervention_effects_results_stratified_arms.RDS"))

res_combined_arm_combined <- bind_rows(res, res_untargeted_metabolomics,
                                       res_microbiome, res_protien)
saveRDS(res_combined_arm_combined,file=paste0(here::here(),"/results/combined_intervention_effects_results_combined_arms.RDS"))

# Native-unit (unscaled) estimates exist only for the targeted panels. The
# untargeted, microbiome and proteomics slots are NULL, so bind_rows() keeps the
# targeted rows alone; they are joined onto the scaled results below.
res_untargeted_metabolomics_unscaled=NULL
res_untargeted_metabolomics_unscaled_arm_strat=NULL
res_microbiome_unscaled=NULL
res_microbiome_unscaled_arm_strat=NULL
res_unscaled_protien=NULL
res_unscaled_protien_arm_strat=NULL
res_combined_arm_combined_unscaled  <- bind_rows(res_unscaled, res_untargeted_metabolomics_unscaled, res_microbiome_unscaled, res_unscaled_protien) %>%
  ungroup() %>%
  select(measure,studytime, contrast, biomarker , est,   cil,   ciu, pval, pval_adj) %>%
  rename(est_unscaled=est, 
         cil_unscaled=cil,   
         ciu_unscaled=ciu, 
         pval_unscaled=pval, 
         pval_adj_unscaled=pval_adj)

res_combined_arm_strat_unscaled <- bind_rows(res_unscaled_arm_strat, res_untargeted_metabolomics_unscaled_arm_strat, res_microbiome_unscaled_arm_strat, res_unscaled_protien_arm_strat) %>%
  ungroup() %>%
  select(measure,studytime, contrast, biomarker , est,   cil,   ciu, pval, pval_adj) %>%
  rename(est_unscaled=est, 
         cil_unscaled=cil,   
         ciu_unscaled=ciu, 
         pval_unscaled=pval, 
         pval_adj_unscaled=pval_adj)

# Percent difference in the arm-stratified native-unit arm means:
# (intervention mean - Control mean) / Control mean * 100.
per_imp_df <- res_combined_arm_strat_unscaled %>% filter(measure!="ATE") %>% group_by(studytime, biomarker) %>%
  mutate(contrast=factor(contrast, levels=c("Control","Nico","BEP"))) %>% arrange(studytime, biomarker, contrast) %>%
  mutate(perc_imp=(last(est_unscaled)-first(est_unscaled))/first(est_unscaled) * 100) %>% 
  filter(contrast!="Control") %>%
  select(studytime, biomarker, contrast, perc_imp)

# Combined-arms ATE subsets by outcome panel: primary_* feed Fig 2,
# secondary_* feed Fig S3.
head(res_combined_arm_combined)
res_combined_arm_combined <- res_combined_arm_combined %>% filter(measure=="ATE")
res_combined_arm_combined <- left_join(res_combined_arm_combined,res_combined_arm_combined_unscaled %>% filter(measure=="ATE"), by=c("studytime", "contrast","biomarker","measure"))
res_combined_arm_combined <- left_join(res_combined_arm_combined,per_imp_df, by=c("studytime", "contrast","biomarker"))

res = res_combined_arm_combined %>% select(study, visit, contrast, category , biomarker, label_f,  est, cil, ciu, pval, pval_adj, sig, sigFDR, est_unscaled, cil_unscaled, ciu_unscaled, pval_adj_unscaled, perc_imp )
res$biomarker = str_to_lower(res$biomarker)

names(all_milk_components)
res_macro = res %>% filter(biomarker %in% c(all_milk_components$macro, "total carbohydrate"))
length(unique(res_macro$biomarker))==length(all_milk_components$macro)

res_micro = res %>% filter(biomarker %in% c(all_milk_components$micro))
length(unique(res_micro$biomarker))==length(all_milk_components$micro)

res_bvit = res %>% filter(biomarker %in% c(all_milk_components$bvit))
length(unique(res_bvit$biomarker))==length(all_milk_components$bvit)

write.csv(res_macro,file=paste0(here::here(),"/results/subsetted results/primary_macro.csv"))
write.csv(res_micro,file=paste0(here::here(),"/results/subsetted results/primary_micro.csv"))
write.csv(res_bvit,file=paste0(here::here(),"/results/subsetted results/primary_bvit.csv"))


head(res_macro)

res_hmo = res %>% filter(biomarker %in% c(all_milk_components$hmo))
length(unique(res_hmo$biomarker))==length(all_milk_components$hmo)

res_bioactives = res %>% filter(biomarker %in% c(all_milk_components$protein))
length(unique(res_bioactives$biomarker))==length(all_milk_components$protein)

write.csv(res_hmo,file=paste0(here::here(),"/results/subsetted results/secondary_hmo.csv"))
write.csv(res_bioactives,file=paste0(here::here(),"/results/subsetted results/secondary_bioactives.csv"))

res_metabolomics = res %>% filter(biomarker %in% c(all_milk_components$metabolomics))
length(unique(res_metabolomics$biomarker))==length(all_milk_components$metabolomics)

write.csv(res_metabolomics,file=paste0(here::here(),"/results/subsetted results/tertiary_targeted_metabolomics.csv"))

# Same subsets from the arm-stratified results.
head(res_combined_arm_strat)
res_combined_arm_strat <- res_combined_arm_strat %>% filter(measure=="ATE")
# ATE rows only on the unscaled side, as in the pooled join above; without the filter each
# effect row was duplicated by the arm-mean (MN) rows of the unscaled results.
res_combined_arm_strat <- left_join(res_combined_arm_strat,res_combined_arm_strat_unscaled %>% filter(measure=="ATE"), by=c("studytime","contrast","biomarker","measure"))

res = res_combined_arm_strat %>% select(study, visit, contrast, category , biomarker, label_f,  est, cil, ciu, pval, pval_adj, sig, sigFDR, est_unscaled, cil_unscaled, ciu_unscaled, pval_adj_unscaled)
res$biomarker = str_to_lower(res$biomarker)

names(all_milk_components)
res_macro = res %>% filter(biomarker %in% c(all_milk_components$macro, "total carbohydrate"))
length(unique(res_macro$biomarker))==length(all_milk_components$macro)

res_micro = res %>% filter(biomarker %in% c(all_milk_components$micro))
length(unique(res_micro$biomarker))==length(all_milk_components$micro)

res_bvit = res %>% filter(biomarker %in% c(all_milk_components$bvit))
length(unique(res_bvit$biomarker))==length(all_milk_components$bvit)

write.csv(res_macro,file=paste0(here::here(),"/results/subsetted results/primary_macro_arm_strat.csv"))
write.csv(res_micro,file=paste0(here::here(),"/results/subsetted results/primary_micro_arm_strat.csv"))
write.csv(res_bvit,file=paste0(here::here(),"/results/subsetted results/primary_bvit_arm_strat.csv"))

head(res_macro)

res_hmo = res %>% filter(biomarker %in% c(all_milk_components$hmo))
length(unique(res_hmo$biomarker))==length(all_milk_components$hmo)

res_bioactives = res %>% filter(biomarker %in% c(all_milk_components$protein))
length(unique(res_bioactives$biomarker))==length(all_milk_components$protein)

write.csv(res_hmo,file=paste0(here::here(),"/results/subsetted results/secondary_hmo_arm_strat.csv"))
write.csv(res_bioactives,file=paste0(here::here(),"/results/subsetted results/secondary_bioactives_arm_strat.csv"))


res_metabolomics = res %>% filter(biomarker %in% c(all_milk_components$metabolomics))
length(unique(res_metabolomics$biomarker))==length(all_milk_components$metabolomics)

write.csv(res_metabolomics,file=paste0(here::here(),"/results/subsetted results/tertiary_targeted_metabolomics_arm_strat.csv"))
