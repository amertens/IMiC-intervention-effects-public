# =============================================================================
# clean_blood_results.R
#
# Tidies the raw bioTMLE output of both script-12 variants (arm-stratified and
# combined-arms) into the long format the milk results use, with Benjamini-Hochberg
# FDR per dataset x visit, as stated in the Methods for the blood analyses. Each
# results object is a named list of per-compartment bioTMLE outputs (MaternalPlasma,
# VamsPrenatal, VamsPostnatalMaternal, VamsPostnatalInfant, ProteomicsDepleted,
# ProteomicsNaive); each run inside it is named "<dataset>-<visit>" (= studytime),
# so grouping by studytime and measure gives FDR per dataset x visit.
# pval_adj_global is pooled across visits within a dataset.
#
# The milk tidier extract_bioTMLE_results() (functions/data_cleaning_functions.R) is
# not used because it relabels milk study times; extract_blood_results()
# (_blood_helpers.R) reuses the same FDR core.
#
# Inputs : results/blood_compartment_[adjusted_]intervention_effects_results.RDS
#          results/blood_compartment_[adjusted_]combined_arms_intervention_effects_results.RDS
# Outputs: the same file names with a _clean suffix (read by combine_blood_results.R
#          and the downstream blood scripts).
# [needs restricted data]
# =============================================================================

# 0-config.R loads extract_res(), ci_to_pvalue() and data.table, which the helper needs.
if (!exists("BLOOD_ORCHESTRATED")) { rm(list = ls()); source(paste0(here::here(), "/src/0-config.R")) }
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
