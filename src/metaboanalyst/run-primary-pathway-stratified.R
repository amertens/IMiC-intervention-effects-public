# run-primary-pathway-stratified.R, arm-STRATIFIED primary via PATHWAY analysis.
#
# Sibling of run-primary-pathway-compare.R (which runs the combined arms). This
# runs the SAME pathway (impact) module on the arm-STRATIFIED contrasts so the
# online supplement's Fig 3B can show the per-arm pathway-impact results, matching
# the arm-stratified Fig 3A volcano grid.
#
#   Combined  (main figure)      -> results/metaboanalyst/primary_combined_pathway/
#   Stratified (this, supplement)-> results/metaboanalyst/primary_stratified_pathway/
#
# One non-directional cell per study x timepoint x arm-contrast (sigFDR==1 features,
# positive and negative analysed together), identical construction to the combined
# run; only the input RDS (stratified arms) and output directory differ.
suppressMessages({ library(dplyr) })
source("src/metaboanalyst/R/run-outcome-group.R")

STRATIFIED_RDS <- "results/combined_intervention_effects_results_stratified_arms.RDS"

run_primary_pathway_stratified <- function(write = TRUE, max_cells = Inf) {
  stratified <- readRDS(STRATIFIED_RDS)
  res <- run_outcome_group(stratified, "primary", module = "pathway",
                           arm_set = "stratified_pathway",
                           query_mode = "combined_sigfdr",
                           write = write, max_cells = max_cells)
  res$results
}

if (sys.nframe() == 0) invisible(run_primary_pathway_stratified(write = TRUE))
