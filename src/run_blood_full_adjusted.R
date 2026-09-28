# =============================================================================
# run_blood_full_adjusted.R
#
# Single-session runner for the covariate-adjusted blood-compartment ATEs
# (MISAME-III), the analysis of record for the blood results (Tables S8-S10,
# Fig. 6D), with BLOOD_ADJUST = TRUE: the milk-pipeline
# covariates (Wvars) are attached through the milk-ID bridge (blood idBiospe = milk
# BMID number -> subjid, 1:1 -> cleaned baseline covariates), and outputs carry an
# adjusted_ tag. About 290 of 309 blood
# samples have a milk/baseline match; the rest are dropped from the adjusted models
# (reported per compartment).
#
# Run from the repo root, then run src/run_blood_adjusted_downstream.R:
#   Rscript src/run_blood_full_adjusted.R
# Outputs: results/blood_compartment_adjusted_intervention_effects_results[_clean].RDS
#          results/blood_compartment_adjusted_combined_arms_intervention_effects_results[_clean].RDS
# [needs restricted data]
# =============================================================================

rm(list = ls())
source(paste0(here::here(), "/src/0-config.R"))

BLOOD_ORCHESTRATED <- TRUE
BLOOD_ADJUST       <- TRUE      # covariate adjustment (FALSE gives the unadjusted companion run)
SETUP_N_FEATURES   <- NULL      # all features (set a number for a quick test run)
N_FEATURES_SUBSET  <- NULL
BLOOD_WORKERS      <- 4         # cap PSOCK workers (31 GB RAM; an unbounded count runs out of memory on the VAMS)
# SL.glm for all compartments: with about 18 covariates, the full milk library on
# the adjusted proteomics took several days.
BLOOD_SL_LIB_PROT  <- c("SL.glm")
BLOOD_SL_LIB_METAB <- c("SL.glm")
here_ <- here::here()

# Protect the config objects and control variables from the between-step cleanup.
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
