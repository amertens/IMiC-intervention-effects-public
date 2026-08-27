# =============================================================================
# run_blood_adjusted_downstream.R
#
# Re-runs the whole DOWNSTREAM blood pipeline on the ADJUSTED results, after
# run_blood_full_adjusted.R has produced the adjusted_ clean files. Everything is
# written to *_adjusted-tagged outputs, so the unadjusted results are preserved.
#
# Steps: combine master -> cross-compartment proteome (13) -> metabolite m/z-RT
# matcher (14) -> mummichog on adjusted blood (15) -> cross-compartment pathway
# figure (17, if present).
#
# Run (after the adjusted run finishes):  Rscript src/run_blood_adjusted_downstream.R
# =============================================================================

root <- paste0(here::here(), "/")
need <- paste0(root, "results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")
if (!file.exists(need))
  stop("Adjusted blood results not found yet:\n  ", need,
       "\nWait for run_blood_full_adjusted.R to finish, then re-run this.")

BLOOD_ADJUST     <- TRUE                                   # -> 13/14/15 read adjusted blood, write _adjusted
RUN_MUMMICHOG    <- TRUE                                   # mummichog on adjusted blood
USE_CONDA        <- TRUE
CONDA_CMD        <- Sys.getenv("IMIC_CONDA_CMD", "C:/Users/andre/miniconda3/Scripts/conda.exe")
CONDA_ENV        <- "mummichog"
MUMMICHOG_ENGINE <- "cli"

src <- paste0(root, "src/2 analysis/")
step <- function(rel) { cat("\n>>>", rel, format(Sys.time()), "\n"); source(paste0(src, rel), local = FALSE) }

step("combine_blood_results.R")                            # master table now includes adjustment=="adjusted"
step("13-cross-compartment-proteome-overlap.R")           # -> *_adjusted.csv/.png
step("14-cross-compartment-metabolite-mzrt-match.R")      # -> *_adjusted.csv/.png
step("15-blood-mummichog-pathway-analysis.R")             # -> mummichog_output_adjusted/, blood_mummichog_pathways_adjusted.*
if (file.exists(paste0(src, "17-cross-compartment-pathway-figure.R")))
  step("17-cross-compartment-pathway-figure.R")           # refresh milk-vs-blood pathway figure on adjusted

cat("\n===== ADJUSTED DOWNSTREAM COMPLETE =====\n")
