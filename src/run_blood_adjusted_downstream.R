# =============================================================================
# run_blood_adjusted_downstream.R
#
# Runs the downstream blood, cross-compartment and annotation scripts on the
# covariate-adjusted blood results, after run_blood_full_adjusted.R has written the
# adjusted_ clean files. Scripts are sourced in dependency order in one session, and
# the global environment is cleared between steps so that a function defined by one
# script (for example annotate() in script 23) cannot mask a package function used by
# a later one. Steps 15, 16 and 18 call Mummichog and need the conda environment
# "mummichog"; set IMIC_CONDA_CMD to the conda executable if conda is not on PATH.
#
# Run from the repo root (after the adjusted run finishes):
#   Rscript src/run_blood_adjusted_downstream.R
# The src/metaboanalyst/ scripts (run-untargeted-msea.R, run-proteomics-go.R,
# run-compartment-pathway.R) and the figure scripts then take these outputs into
# Fig. 6 and Tables S5, S7 and S11.
# [needs restricted data]
# =============================================================================

root <- paste0(here::here(), "/")
need <- paste0(root, "results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")
if (!file.exists(need))
  stop("Adjusted blood results not found yet:\n  ", need,
       "\nWait for run_blood_full_adjusted.R to finish, then re-run this.")

BLOOD_ADJUST     <- TRUE                                   # 13 and 15 read adjusted blood, write _adjusted outputs
RUN_MUMMICHOG    <- TRUE                                   # 15 runs Mummichog, not only its inputs
USE_CONDA        <- TRUE
CONDA_CMD        <- Sys.getenv("IMIC_CONDA_CMD", "conda")
CONDA_ENV        <- "mummichog"
MUMMICHOG_ENGINE <- "cli"

src <- paste0(root, "src/2 analysis/")
PROTECT <- c(ls(), "PROTECT", "step")                      # settings kept across steps
step <- function(rel) {
  t <- Sys.time()
  cat("\n>>>", rel, format(t), "\n")
  source(paste0(src, rel), local = FALSE)
  cat(sprintf(">>> done in %.1f min\n", as.numeric(difftime(Sys.time(), t, units = "mins"))))
  rm(list = setdiff(ls(envir = .GlobalEnv), PROTECT), envir = .GlobalEnv); invisible(gc())
}

# Blood master table and cross-compartment proteome
step("combine_blood_results.R")                            # -> blood_compartment_all_*.{RDS,csv}
step("13-cross-compartment-proteome-overlap.R")            # -> cross_compartment_proteome_overlap_adjusted.csv
# Mummichog runs (conda)
step("15-blood-mummichog-pathway-analysis.R")              # -> mummichog_output_adjusted/
step("16-milk-mummichog-for-comparison.R")                 # -> mummichog_output/Milk_*
step("18-directional-mummichog.R")                         # -> mummichog_output_directional/
# Direction agreement across compartments
step("20-cross-compartment-threshold-free-panel.R")        # -> cross_compartment_threshold_free_panel.csv
step("29-figure-overall-shift-agreement.R")                # -> 16-of-18 maternal->infant agreement (console)
# Annotation and pathway direction
step("23-annotate-fdr-features.R")                         # needs 15 -> fdr_sig_putative_annotation.csv
step("32-signed-pathway-direction.R")                      # needs 15 -> signed_pathway_direction.csv
step("33-signed-pathway-direction-milk.R")                 # needs 16 -> signed_pathway_direction_milk.csv
step("41-fat-synthesis-timepoint-table.R")                 # needs 33 -> fat_synthesis_timepoint_table.csv
step("47-annotate-milk-features.R")                        # needs 18 -> milk_{nominal,fdr_sig}_putative_annotation.csv
# Supplementary tables and supplement matching
step("46-build-table-s8-cross-compartment.R")              # needs 13 -> Table S8
step("52-blood-class-enrichment-direction-sensitivity.R")  # -> Table S9
step("57-table-s10-temporal-persistence.R")                # needs combine -> Table S10
step("50-blood-mummichog-input.R")                         # -> blood_mummichog_input.csv
step("54-supplement-detection-status.R")                   # needs 50, 23, 47 -> supplement_status_fdr_features.csv
# Milk proteome GO enrichment (independent of the blood results)
step("55-proteomics-go-uniprot.R")                         # -> proteomics_go_uniprot.csv (Fig. 6C, Table S7)

cat("\n===== ADJUSTED DOWNSTREAM COMPLETE =====\n")
