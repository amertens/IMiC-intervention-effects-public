# run-primary-pathway-compare.R — combined-arm primary via PATHWAY analysis.
#
# Why this exists: the manuscript's automated primary enrichment (run-primary.R)
# uses the Over-Representation Analysis (ORA) module. Trenton's ORIGINAL published
# primary figure instead used the Pathway (impact) analysis module on the combined
# arms. To compare the two on the same footing, this driver runs the combined-arm
# primary cells through the pathway module (pathora) and writes them to their own
# directory so they never overwrite the ORA results.
#
#   ORA (main analysis)      -> results/metaboanalyst/primary_combined/
#   Pathway (this, compare)  -> results/metaboanalyst/primary_combined_pathway/
#
# The two runs use the SAME per-cell query feature list; only the enrichment
# module differs (msetora vs pathora), which is the apples-to-apples comparison.
# Pathway analysis takes no reference metabolome (that is an ORA-only setting).
suppressMessages({ library(dplyr) })
source("src/metaboanalyst/R/run-outcome-group.R")

# Repo's own ATE output (byte-identical to Trenton's), for end-to-end raw -> ATE -> MetaboAnalyst.
COMBINED_RDS <- "results/combined_intervention_effects_results_combined_arms.RDS"

run_primary_pathway <- function(write = TRUE, max_cells = Inf) {
  combined <- readRDS(COMBINED_RDS)
  # arm_set label "combined_pathway" only sets the output directory / labels;
  # it does not change the analysis (all studies, combined-arm contrasts).
  # Fig 3B construction (matches Trenton's "Primary Outcomes (Pathway Analysis)"):
  # one non-directional cell per study x timepoint, query = sigFDR==1 features with
  # positive and negative estimates analysed together (query_mode="combined_sigfdr").
  res <- run_outcome_group(combined, "primary", module = "pathway",
                           arm_set = "combined_pathway",
                           query_mode = "combined_sigfdr",
                           write = write, max_cells = max_cells)
  res$results
}

if (sys.nframe() == 0) invisible(run_primary_pathway(write = TRUE))
