# =============================================================================
# run_blood_adjusted_recover.R
#
# The adjusted run (run_blood_full_adjusted.R) finished its STRATIFIED step and
# saved blood_compartment_adjusted_intervention_effects_results.RDS, then OOM-died
# at the start of the COMBINED step (4 PSOCK workers x 38k-feature VAMS x 18
# covariates blew past 31 GB). This recovers from there:
#   combined (2 workers) -> clean -> combine master -> cross-compartment (13,14) ->
#   mummichog (15) -> pathway figure (17), all ADJUSTED / *_adjusted-tagged.
#
# Uses the on-disk adjusted merged_blood_datasets.RDS from the 08:25 prep (verified
# to carry the Wvars covariates) -- no need to redo prep or stratified.
#
# Run:  Rscript src/run_blood_adjusted_recover.R   (long; ~5-6h, mostly combined)
# =============================================================================

rm(list = ls()); source(paste0(here::here(), "/src/0-config.R"))

BLOOD_ORCHESTRATED <- TRUE
BLOOD_ADJUST       <- TRUE                       # read adjusted blood, write *_adjusted
BLOOD_WORKERS      <- 2                           # <-- fewer workers: avoid the combined-step OOM
N_FEATURES_SUBSET  <- NULL                        # full features
BLOOD_SL_LIB_METAB <- c("SL.glm")
BLOOD_SL_LIB_PROT  <- c("SL.glm")
RUN_MUMMICHOG      <- TRUE
USE_CONDA          <- TRUE
CONDA_CMD          <- Sys.getenv("IMIC_CONDA_CMD", "C:/Users/andre/miniconda3/Scripts/conda.exe")
CONDA_ENV          <- "mummichog"
MUMMICHOG_ENGINE   <- "cli"
here_ <- here::here()

step <- function(rel) {
  t <- Sys.time(); cat("\n>>>", rel, "start", format(t), "\n")
  source(paste0(here_, "/", rel), local = FALSE)
  cat(">>>", rel, sprintf("done %.1f min\n", as.numeric(difftime(Sys.time(), t, units = "mins"))))
  gc()
}

step("src/2 analysis/12-blood-compartment-intervention-effects_combined_arms.R")  # adjusted combined raw
step("src/2 analysis/clean_blood_results.R")                                       # adjusted clean (strat+comb)
step("src/2 analysis/combine_blood_results.R")                                     # master incl. adjusted rows
step("src/2 analysis/13-cross-compartment-proteome-overlap.R")                    # *_adjusted
step("src/2 analysis/14-cross-compartment-metabolite-mzrt-match.R")               # *_adjusted
step("src/2 analysis/15-blood-mummichog-pathway-analysis.R")                      # mummichog on adjusted
step("src/3 visualizations/17-cross-compartment-pathway-figure.R")               # adjusted pathway figure

cat("\n===== ADJUSTED RECOVERY + DOWNSTREAM COMPLETE =====\n")
