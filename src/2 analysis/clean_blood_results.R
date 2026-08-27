
# =============================================================================
# clean_blood_results.R
#
# Tidies raw bioTMLE output from the blood-compartment ATE scripts (stratified +
# combined-arms) into the long format the milk results use.
#
# FDR: BH GROUPED BY VISIT x DATASET (x measure). Each run is named
# "<DatasetLabel>-<visit>" (= studytime), so group_by(studytime, measure) ==
# per dataset x visit. pval_adj_global is pooled across visits within a dataset.
#
# Results objects are now a NAMED LIST of per-compartment bioTMLE outputs
# (MaternalPlasma, VamsPrenatal, VamsPostnatalMaternal, VamsPostnatalInfant,
# ProteomicsDepleted, ProteomicsNaive). We do NOT use extract_bioTMLE_results()
# (it hardcodes milk studytime relabeling); extract_blood_results() reuses the
# same FDR core.
# =============================================================================

if (!exists("BLOOD_ORCHESTRATED")) { rm(list = ls()); source(paste0(here::here(), "/src/0-config.R")) }
# 0-config.R provides extract_res, ci_to_pvalue, data.table.
# extract_blood_results() (raw bioTMLE list -> long, FDR per visit x dataset) is the
# canonical version in the shared helpers.
source(paste0(here::here(), "/src/2 analysis/_blood_helpers.R"))

clean_blood_set <- function(results) {
  dplyr::bind_rows(Map(extract_blood_results, results, names(results)))
}

if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- FALSE
.tag <- if (BLOOD_ADJUST) "adjusted_" else ""   # match the analysis output tag

# --- arm-stratified ----------------------------------------------------------
res_strat <- readRDS(paste0(here::here(),
  "/results/blood_compartment_", .tag, "intervention_effects_results.RDS"))
saveRDS(clean_blood_set(res_strat), file = paste0(here::here(),
  "/results/blood_compartment_", .tag, "intervention_effects_results_clean.RDS"))

# --- combined arms -----------------------------------------------------------
res_comb <- readRDS(paste0(here::here(),
  "/results/blood_compartment_", .tag, "combined_arms_intervention_effects_results.RDS"))
saveRDS(clean_blood_set(res_comb), file = paste0(here::here(),
  "/results/blood_compartment_", .tag, "combined_arms_intervention_effects_results_clean.RDS"))

# Next: join cleaned blood ATEs to cleaned milk results on aligned features
# (IMiC_alignment.csv for metabolomics; shared UniProt for proteomics) and
# compare ATE direction/significance across milk / maternal blood / infant blood.
