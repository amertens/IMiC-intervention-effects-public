# =============================================================================
# run_blood_full_adjusted.R  --  single-session ADJUSTED blood-compartment runner
#
# Same as run_blood_full.R but with covariate adjustment ON (BLOOD_ADJUST = TRUE):
# attaches the milk-pipeline Wvars via the milk-ID bridge
#   blood idBiospe == milk BMID number -> subjid (1:1) -> clean_baseline_covariates
# and writes to ADJUSTED-tagged outputs, so it does NOT overwrite the unadjusted
# results:
#   results/blood_compartment_adjusted_intervention_effects_results[_clean].RDS
#   results/blood_compartment_adjusted_combined_arms_intervention_effects_results[_clean].RDS
#
# Run:  Rscript src/run_blood_full_adjusted.R
# Do NOT run concurrently with run_blood_full.R (both rebuild merged_blood_datasets.RDS
# and compete for memory) -- launch this AFTER the unadjusted run finishes.
#
# Coverage note: ~290/309 blood samples carry a milk/baseline match; the ~19 without
# a milk sample are dropped from the adjusted models (reported per compartment).
# =============================================================================

rm(list = ls())
source(paste0(here::here(), "/src/0-config.R"))

BLOOD_ORCHESTRATED <- TRUE
BLOOD_ADJUST       <- TRUE      # <-- the only difference from run_blood_full.R
SETUP_N_FEATURES   <- NULL      # full features; set to 100 for a quick adjusted setup pass
N_FEATURES_SUBSET  <- NULL
BLOOD_WORKERS      <- 4         # cap PSOCK workers (31GB RAM; unbounded OOMs on 38k-feat VAMS)
# SL.glm for ALL compartments so the adjusted run finishes in ~12-24h (one uptime
# window). The full library on adjusted proteomics (~18 covariates) proved too slow
# / multi-day; revisit it on HPC. (Override BLOOD_SL_LIB_PROT to re-enable.)
BLOOD_SL_LIB_PROT  <- c("SL.glm")
BLOOD_SL_LIB_METAB <- c("SL.glm")
here_ <- here::here()

PROTECT <- c(ls(), "PROTECT", "step", "t0")

step <- function(label, relpath) {
  t <- Sys.time()
  cat(sprintf("\n>>> [%s] start %s\n", label, format(t)))
  source(paste0(here_, relpath), local = FALSE)
  cat(sprintf(">>> [%s] done in %.1f min\n", label,
              as.numeric(difftime(Sys.time(), t, units = "mins"))))
  rm(list = setdiff(ls(envir = .GlobalEnv), PROTECT), envir = .GlobalEnv); gc()
}

t0 <- Sys.time()
step("prep(adj)",       "/src/1 data prep/7-blood-compartment-prep.R")
step("stratified(adj)", "/src/2 analysis/12-blood-compartment-intervention-effects.R")
step("combined(adj)",   "/src/2 analysis/12-blood-compartment-intervention-effects_combined_arms.R")
step("clean(adj)",      "/src/2 analysis/clean_blood_results.R")
cat(sprintf("\n===== ALL_BLOOD_ADJUSTED_DONE in %.1f min =====\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
