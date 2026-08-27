# =============================================================================
# 24-fdr-updown-feature-lists.R
#
# Trenton's request: a clean export of the up- and down-regulated FDR-significant
# metabolite features, per dataset x timepoint x treatment contrast, with mass,
# retention time and (where available) a name. This is the raw material the
# matching in 19b consumes and the small lists we hand to Kim.
#
# Out: results/fdr_sig_updown_lists.csv
# =============================================================================

suppressMessages({library(data.table)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))

# m/z-RT + names (postnatal VAMS named catalogue; milk pipeline description; rLC unnamed)
vm  <- fread(paste0(root, "data/additional datasets/metabolite_description_vam_with_global_id.csv"))
vmz <- data.table(feature = toupper(vm$feature_label), mz = vm$mz, rt = vm$rt_minute,
                  name = ifelse(is.na(vm$ID) | vm$ID == "", NA_character_, as.character(vm$ID)))
plz <- mzrt_rlc("ProcessedDataMISAME3_plasma.csv"); pvz <- mzrt_rlc("ProcessedDataMISAME3_VAMS.csv"); mkz <- mzrt_milk(misame_only = FALSE)
mzlk <- unique(rbindlist(list(vmz[, .(feature, mz, rt)], plz[, .(feature, mz, rt)],
                              pvz[, .(feature, mz, rt)], mkz[, .(feature, mz, rt)]), use.names = TRUE))
mz_of <- function(f) mzlk$mz[match(toupper(f), mzlk$feature)]
rt_of <- function(f) mzlk$rt[match(toupper(f), mzlk$feature)]
nm_of <- function(f) vmz$name[match(toupper(f), vmz$feature)]

milk  <- as.data.table(readRDS(paste0(root,"results/adjusted_intervention_effects_res_untargeted_metabolomics_clean.RDS")))[measure=="ATE" & study=="Misame"]
blood <- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_intervention_effects_results_clean.RDS")))[measure=="ATE"]
milkC <- as.data.table(readRDS(paste0(root,"results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")))[measure=="ATE" & study=="Misame"]
bloodC<- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))[measure=="ATE"]

# one tidy FDR-sig table from any clean results, with a dataset label and a `contrast` column
tidy <- function(d, dataset_label) {
  d <- copy(d); d[, feature := toupper(biomarker)]
  d <- d[sigFDR == 1][order(pval)]
  if (!nrow(d)) return(NULL)
  dt <- if ("dataset" %in% names(d)) d$dataset else dataset_label
  data.table(dataset = dt, visit = d$visit, contrast = d$contrast,
             direction = ifelse(d$est > 0, "up", "down"), feature = d$feature,
             est = round(d$est, 2), q = signif(d$pval_adj, 2),
             mz = round(mz_of(d$feature), 4), rt = round(rt_of(d$feature), 3), name = nm_of(d$feature))
}
milkC[, contrast := "any-postnatal-BEP"]; bloodC[, contrast := "any-postnatal-BEP"]
res <- rbindlist(list(
  tidy(milk,   "Milk (untargeted)"),          # milk, arm-stratified contrasts
  tidy(milkC,  "Milk (untargeted)"),          # milk, any-postnatal-BEP
  tidy(blood,  NA),                            # blood, arm-stratified (dataset col present)
  tidy(bloodC, NA)                             # blood, any-postnatal-BEP
), use.names = TRUE, fill = TRUE)
fwrite(res, paste0(root, "results/fdr_sig_updown_lists.csv"))
cat("Saved results/fdr_sig_updown_lists.csv :", nrow(res), "FDR-sig feature rows\n")
cat("by dataset x direction:\n"); print(res[, .N, by = .(dataset, direction)][order(dataset)])
cat("named features:", sum(!is.na(res$name)), "\n")
