# =============================================================================
# src/2 analysis/clean_results.R
#
# Reads:  metadata/milk_component.Rdata
#         results/adjusted_combined_arms_HMtraj_intervention_effects_results.RDS
#         results/adjusted_combined_arms_intervention_effects_results.RDS
#         results/adjusted_combined_arms_intervention_effects_results_unscaled.RDS
#         results/adjusted_combined_arms_tertiary_intervention_effects_results_unscaled.RDS
#         results/adjusted_HMtraj_intervention_effects_results.RDS
#         results/adjusted_intervention_effects_results.RDS
#         results/adjusted_intervention_effects_results_unscaled.RDS
#         results/microbiome_intervention_effects_results.RDS
#         results/microbiome_intervention_effects_results_arm_strat.RDS
#         results/proteomics_intervention_effects_results.RDS
#         results/proteomics_intervention_effects_results_combined_arms.RDS
#         results/temp_vital_adjusted_combined_arms_tertiary_intervention_effects_results.RDS
#         results/unadjusted_intervention_effects_results.RDS
# Writes: results/adjusted_combined_arms_intervention_effects_proteomics_results_clean.RDS
#         results/adjusted_combined_arms_intervention_effects_results_clean.csv
#         results/adjusted_combined_arms_intervention_effects_results_clean.RDS
#         results/adjusted_combined_arms_intervention_effects_traj_results_clean.RDS
#         results/adjusted_combined_arms_intervention_effects_unscaled_results_clean.RDS
#         results/adjusted_combined_arms_intervention_effects_untargeted_results_clean.RDS
#         results/adjusted_intervention_effects_res_untargeted_metabolomics_clean.RDS
#         results/adjusted_intervention_effects_results_clean.RDS
#         results/adjusted_intervention_effects_traj_results_clean.RDS
#         results/adjusted_intervention_effects_unscaled_results_clean.RDS
#         results/combined_intervention_effects_results_combined_arms.RDS
#         results/combined_intervention_effects_results_combined_arms_ATE_sharing.csv
#         results/combined_intervention_effects_results_combined_arms_ATE_sharing.RDS
#         results/combined_intervention_effects_results_combined_arms_ATE_sharing_trenton.csv
#         results/combined_intervention_effects_results_combined_arms_ATE_sharing_trenton.RDS
#         results/combined_intervention_effects_results_stratified_arms.RDS
#         results/combined_intervention_effects_results_stratified_arms_ATE_sharing.csv
#         results/combined_intervention_effects_results_stratified_arms_ATE_sharing.RDS
#         results/combined_intervention_effects_results_stratified_arms_ATE_sharing_trenton.csv
#         results/combined_intervention_effects_results_stratified_arms_ATE_sharing_trenton.RDS
#         results/results_subset.csv
#         results/subsetted results/primary_bvit.csv
#         results/subsetted results/primary_bvit_arm_strat.csv
#         results/subsetted results/primary_macro.csv
#         results/subsetted results/primary_macro_arm_strat.csv
#         results/subsetted results/primary_micro.csv
#         results/subsetted results/primary_micro_arm_strat.csv
#         results/subsetted results/secondary_bioactives.csv
#         results/subsetted results/secondary_bioactives_arm_strat.csv
#         results/subsetted results/secondary_hmo.csv
#         results/subsetted results/secondary_hmo_arm_strat.csv
#         results/subsetted results/tertiary_targeted_metabolomics.csv
#         results/subsetted results/tertiary_targeted_metabolomics_arm_strat.csv
#         results/unadjusted_intervention_effects_results_clean.RDS
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



res <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results.RDS"))
res_vital_tertiary <- readRDS(file=paste0(here::here(),"/results/temp_vital_adjusted_combined_arms_tertiary_intervention_effects_results.RDS"))
res$res_tertiary$res$`Vital-40`$res <- res_vital_tertiary$res[[1]]$res
res$res_tertiary$res$`Vital-56`$res <- res_vital_tertiary$res[[2]]$res

res_unscaled <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_tertiary_intervention_effects_results_unscaled.RDS"))
res_traj <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_HMtraj_intervention_effects_results.RDS"))

res_untargeted_metabolomics <- readRDS(paste0(here::here(),"/large-file-results/metabalomics_intervention_effects_results.RDS"))



#res_unadjusted <- readRDS(file=paste0(here::here(),"/results/unadjusted_combined_arms_intervention_effects_results.RDS"))
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
#res_unadjusted <- extract_bioTMLE_results(res_unadjusted)
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
head(res_protien)

#Check B status between misame and vital
temp=res_unscaled %>% filter(study  !="Elicit", measure =="MN", contrast=="Control", outcome_group == "primary", grepl("B", biomarker)|biomarker=="Nufa"|biomarker=="Nam"|biomarker=="Ribo"|biomarker=="Pa") 
dim(temp)
table(temp$biomarker)
temp %>% group_by(biomarker, study, visit) %>% summarize(mean_est=mean(est)) %>% as.data.frame()
ggplot(temp, aes(x=visit , y=est, color=study)) + geom_point() + geom_line(aes(group=study)) + facet_wrap(~biomarker, scales="free") + theme_bw()

temp=res_unscaled %>% filter(study  !="Elicit", measure !="MN", outcome_group == "primary", grepl("B", biomarker)|biomarker=="Nufa"|biomarker=="Nam"|biomarker=="Ribo"|biomarker=="Pa") 
dim(temp)
table(temp$biomarker)
temp %>% group_by(biomarker, study, visit) %>% summarize(mean_est=mean(est)) %>% as.data.frame()
ggplot(temp, aes(x=visit , y=est, color=study)) + geom_point() + geom_line(aes(group=study)) + facet_wrap(~biomarker, scales="free") + theme_bw()


temp=res %>% filter(study  !="Elicit", measure !="MN", outcome_group == "primary", grepl("B", biomarker)|biomarker=="Nufa"|biomarker=="Nam"|biomarker=="Ribo"|biomarker=="Pa") 
dim(temp)
table(temp$biomarker)
temp %>% group_by(biomarker, study, visit) %>% summarize(mean_est=mean(est)) %>% as.data.frame()
ggplot(temp, aes(x=visit , y=est, color=study)) + geom_point() + geom_line(aes(group=study)) + facet_wrap(~biomarker, scales="free") + theme_bw()


#res_traj <- extract_bioTMLE_results(res_traj)
#res_untargeted_metabolomics <- extract_bioTMLE_results(res_untargeted_metabolomics, single_group = TRUE)

saveRDS(res,file=paste0(here::here(), "/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))
#saveRDS(res_unadjusted,file=paste0(here::here(), "/results/unadjusted_combined_arms_intervention_effects_results_clean.RDS"))
saveRDS(res_arm_strat,file=paste0(here::here(), "/results/adjusted_intervention_effects_results_clean.RDS"))
saveRDS(res_unadjusted_arm_strat,file=paste0(here::here(), "/results/unadjusted_intervention_effects_results_clean.RDS"))
saveRDS(res_traj_arm_strat,file=paste0(here::here(), "/results/adjusted_combined_arms_intervention_effects_traj_results_clean.RDS"))
saveRDS(res_traj,file=paste0(here::here(), "/results/adjusted_intervention_effects_traj_results_clean.RDS"))
saveRDS(res_unscaled,file=paste0(here::here(), "/results/adjusted_combined_arms_intervention_effects_unscaled_results_clean.RDS"))
saveRDS(res_unscaled_arm_strat,file=paste0(here::here(), "/results/adjusted_intervention_effects_unscaled_results_clean.RDS"))

saveRDS(res_protien,file=paste0(here::here(), "/results/adjusted_combined_arms_intervention_effects_proteomics_results_clean.RDS"))
saveRDS(res_untargeted_metabolomics,file=paste0(here::here(), "/results/adjusted_combined_arms_intervention_effects_untargeted_results_clean.RDS"))



# FIX: save the genuinely *arm-stratified* untargeted object here. Previously this
# line re-saved the *combined-arms* object (`res_untargeted_metabolomics`), so the
# arm-strat filename was byte-identical to the combined-arms file above. The sole
# consumer, 7-metabolomics_plots.R, filters measure=="ATE" and then collapses to
# the max-|effect| row per (studytime, biomarker), a dedup that is only
# meaningful for arm-stratified data, so it was written for this object and is
# corrected (not broken) by this change.
saveRDS(res_untargeted_metabolomics_arm_strat,file=paste0(here::here(), "/results/adjusted_intervention_effects_res_untargeted_metabolomics_clean.RDS"))

## ---- separate ATE (intervention-effect) vs MN (arm-mean) rows --------------
## The high-dimensional clean frames above interleave each ATE estimate with two
## MN (arm-mean) rows (one Control mean + one intervention mean). Any
## "N of M estimates significant" count MUST be computed on ATE rows only:
## counting significance over all rows inflates BOTH numerator and denominator
## (e.g. the untargeted 2,098/1,060,968 and proteome 503/9,003 in-text counts
## were all-row counts; the ATE-only values are 1,347/353,656 and 10/3,001).
## Emit explicit *_ATE.RDS / *_MN.RDS companions so downstream count code never
## has to remember to filter. Existing combined files are retained unchanged for
## backward compatibility (all figure scripts already filter measure=="ATE").
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


saveRDS(res,file=paste0(here::here(),
                        "/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))
write.csv(res,file=paste0(here::here(),
                          "/results/adjusted_combined_arms_intervention_effects_results_clean.csv"))

#subset for results comparison
head(res)
res_sub = res %>% filter()
write.csv(res_sub,file=paste0(here::here(),
                              "/results/results_subset.csv"))


# save combined results 
res_combined_arm_strat <- bind_rows(res_arm_strat, res_untargeted_metabolomics_arm_strat, res_microbiome_arm_strat, res_protien_arm_strat)
saveRDS(res_combined_arm_strat,file=paste0(here::here(),"/results/combined_intervention_effects_results_stratified_arms.RDS"))

res_combined_arm_combined <- bind_rows(res, res_untargeted_metabolomics,
                                       res_microbiome, res_protien)
saveRDS(res_combined_arm_combined,file=paste0(here::here(),"/results/combined_intervention_effects_results_combined_arms.RDS"))


## ---- ATE-only, sharing, modality-tagged merged datasets ---------------------
## A leaner companion to the combined frames above, requested for downstream use:
##   * ATE rows only (MN / arm-mean rows dropped)
##   * uncorrected p-values only (pval / sig kept; all FDR and chi-square columns
##     -- pval_adj, pval_adj_global, sigFDR, chi_pval, chi_pval_adj, chi_sig,
##     chi_sigFDR -- dropped)
##   * a `modality` column naming which assay each combined-modality result came from
## Built in both arm-stratified and combined-arms (merged) versions.
build_ate_sharing <- function(...) {
  drop_cols <- c("pval_adj", "pval_adj_global", "sigFDR",
                 "chi_pval", "chi_pval_adj", "chi_sig", "chi_sigFDR")
  bind_rows(...) %>%
    filter(measure == "ATE") %>%
    select(-any_of(drop_cols))
}

# combined-arms (merged across intervention arms)
res_ate_sharing_combined_arms <- build_ate_sharing(
  res                          %>% mutate(modality = "Targeted milk components"),
  res_untargeted_metabolomics  %>% mutate(modality = "Untargeted metabolomics"),
  res_microbiome               %>% mutate(modality = "Microbiome"),
  res_protien                  %>% mutate(modality = "Proteomics")
)
saveRDS(res_ate_sharing_combined_arms,
        file=paste0(here::here(),"/results/combined_intervention_effects_results_combined_arms_ATE_sharing.RDS"))
write.csv(res_ate_sharing_combined_arms,
          file=paste0(here::here(),"/results/combined_intervention_effects_results_combined_arms_ATE_sharing.csv"), row.names=FALSE)

## ---- Trenton export: minimal collaborator dataset --------------------------
## A trimmed export requested for sharing: ATE point estimate + 95% CI, the raw
## (uncorrected) p-value, the milk component (raw code AND readable label), and
## its study / visit / intervention arm (contrast) / dataset (modality). It
## deliberately drops every FDR column (pval_adj, pval_adj_global, sigFDR), every
## chi-square column (chi_*), and the raw `sig` flag -- raw p-values only.
## Derived from the combined-arms ATE-sharing frame above (already ATE-only,
## uncorrected-p-only, modality-tagged), so it is byte-consistent with it.
trenton_cols <- c("biomarker", "label_f", "study", "visit",
                  "contrast", "modality", "est", "cil", "ciu", "pval")
res_trenton <- res_ate_sharing_combined_arms[, trenton_cols, drop = FALSE]
saveRDS(res_trenton,
        file=paste0(here::here(),"/results/combined_intervention_effects_results_combined_arms_ATE_sharing_trenton.RDS"))
write.csv(res_trenton,
          file=paste0(here::here(),"/results/combined_intervention_effects_results_combined_arms_ATE_sharing_trenton.csv"), row.names=FALSE)
cat(sprintf("[Trenton export] %d rows, cols: %s\n", nrow(res_trenton), paste(trenton_cols, collapse=", ")))

# arm-stratified
res_ate_sharing_arm_strat <- build_ate_sharing(
  res_arm_strat                          %>% mutate(modality = "Targeted milk components"),
  res_untargeted_metabolomics_arm_strat  %>% mutate(modality = "Untargeted metabolomics"),
  res_microbiome_arm_strat               %>% mutate(modality = "Microbiome"),
  res_protien_arm_strat                  %>% mutate(modality = "Proteomics")
)
saveRDS(res_ate_sharing_arm_strat,
        file=paste0(here::here(),"/results/combined_intervention_effects_results_stratified_arms_ATE_sharing.RDS"))
write.csv(res_ate_sharing_arm_strat,
          file=paste0(here::here(),"/results/combined_intervention_effects_results_stratified_arms_ATE_sharing.csv"), row.names=FALSE)

## Trenton export: arm-stratified twin of the trimmed dataset above (same 10 cols).
res_trenton_strat <- res_ate_sharing_arm_strat[, trenton_cols, drop = FALSE]
saveRDS(res_trenton_strat,
        file=paste0(here::here(),"/results/combined_intervention_effects_results_stratified_arms_ATE_sharing_trenton.RDS"))
write.csv(res_trenton_strat,
          file=paste0(here::here(),"/results/combined_intervention_effects_results_stratified_arms_ATE_sharing_trenton.csv"), row.names=FALSE)
cat(sprintf("[Trenton export] arm-strat: %d rows, cols: %s\n", nrow(res_trenton_strat), paste(trenton_cols, collapse=", ")))

cat(sprintf("[ATE sharing] combined-arms: %d rows, %d modalities (%s)\n",
            nrow(res_ate_sharing_combined_arms),
            length(unique(res_ate_sharing_combined_arms$modality)),
            paste(unique(res_ate_sharing_combined_arms$modality), collapse=", ")))
cat(sprintf("[ATE sharing] arm-strat:    %d rows, %d modalities (%s)\n",
            nrow(res_ate_sharing_arm_strat),
            length(unique(res_ate_sharing_arm_strat$modality)),
            paste(unique(res_ate_sharing_arm_strat$modality), collapse=", ")))
cat(sprintf("[ATE sharing] measure levels kept -> combined: %s | arm-strat: %s\n",
            paste(unique(res_ate_sharing_combined_arms$measure), collapse="/"),
            paste(unique(res_ate_sharing_arm_strat$measure), collapse="/")))

res_untargeted_metabolomics_unscaled=NULL #need to add
res_untargeted_metabolomics_unscaled_arm_strat=NULL #need to add
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




#get relative improvements
per_imp_df <- res_combined_arm_strat_unscaled %>% filter(measure!="ATE") %>% group_by(studytime, biomarker) %>%
  mutate(contrast=factor(contrast, levels=c("Control","Nico","BEP"))) %>% arrange(studytime, biomarker, contrast) %>%
  mutate(perc_imp=(last(est_unscaled)-first(est_unscaled))/first(est_unscaled) * 100) %>% 
  filter(contrast!="Control") %>%
  select(studytime, biomarker, contrast, perc_imp)





#save subsets for paper
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



#Arm strat

#save subsets for paper
head(res_combined_arm_strat)
res_combined_arm_strat <- res_combined_arm_strat %>% filter(measure=="ATE")
res_combined_arm_strat <- left_join(res_combined_arm_strat,res_combined_arm_strat_unscaled, by=c("studytime","contrast","biomarker"))

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






