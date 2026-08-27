# run-primary.R — reproduce Trenton's primary-outcome enrichment analysis.
# Thin wrapper over the general run_outcome_group() driver: primary uses the ORA
# module with the curated matched-name reference metabolome (reproduces the
# nicotinate total 7 / hits 5).
#
# Per the 2026-07-30 Andrew x Trenton decision, produce BOTH arm framings for
# every study: combined-arm ORA feeds the MAIN text, stratified-arm ORA feeds
# the SUPPLEMENT. (Previously this was a hybrid -- Elicit/Vital combined,
# Misame stratified -- matching Trenton's original Rmd.) The combined run reads
# every study from the combined-arms RDS; the stratified run reads every study
# from the stratified-arms RDS. Writes results/metaboanalyst/primary_combined/
# and results/metaboanalyst/primary_stratified/.
suppressMessages({ library(dplyr) })
source("src/metaboanalyst/R/config-primary.R")
source("src/metaboanalyst/R/run-outcome-group.R")

COMBINED_RDS   <- "results/combined_intervention_effects_results_combined_arms.RDS"
STRATIFIED_RDS <- "results/combined_intervention_effects_results_stratified_arms.RDS"

# The curated reference metabolome is REQUIRED: it holds MetaboAnalyst-matched
# names that reproduce Trenton's background (nicotinate 7/5). Running without it
# silently changes every ORA p-value, so fail loudly rather than degrade.
.primary_reference <- function() {
  ref_path <- PRIMARY_CONFIG$reference_path
  if (is.null(ref_path) || !file.exists(ref_path)) {
    stop("Primary ORA reference metabolome not found at '", ref_path, "'. ",
         "This curated file is required to reproduce Trenton's background ",
         "(nicotinate total 7 / hits 5); running without it silently changes ",
         "every p-value. Restore the file or fix PRIMARY_CONFIG$reference_path.",
         call. = FALSE)
  }
  trimws(readLines(ref_path, warn = FALSE))
}

run_primary <- function(write = TRUE, max_cells = Inf) {
  ref <- .primary_reference()
  combined   <- readRDS(COMBINED_RDS)
  stratified <- readRDS(STRATIFIED_RDS)

  # All studies, both framings: combined arms -> main text, stratified -> supplement.
  combined_run   <- run_outcome_group(combined,   "primary", module = "ora", arm_set = "combined",
                                      reference_names = ref, write = write, max_cells = max_cells)
  stratified_run <- run_outcome_group(stratified, "primary", module = "ora", arm_set = "stratified",
                                      reference_names = ref, write = write, max_cells = max_cells)

  bind_rows(combined_run$results, stratified_run$results)
}

if (sys.nframe() == 0) invisible(run_primary(write = TRUE))
