# =============================================================================
# run_blood_full.R  --  single-session FULL-feature blood-compartment runner
#
# Loads 0-config.R ONCE, then runs prep -> stratified -> combined -> clean in the
# same R session with the full feature set (overrides the per-script 100-feature
# setup defaults). Frees memory between steps. Prints per-step wall-clock timing.
#
# Run from a shell:
#   Rscript src/run_blood_full.R
# (long job -- the ~38k-feature postnatal VAMS dominates; run detached.)
# Standalone scripts still work unchanged: each only re-loads config / applies the
# 100-feature default when BLOOD_ORCHESTRATED / the override vars are absent.
# =============================================================================

rm(list = ls())
source(paste0(here::here(), "/src/0-config.R"))

BLOOD_ORCHESTRATED <- TRUE     # tells sourced scripts to skip their own rm()+config load
SETUP_N_FEATURES   <- NULL     # prep: read ALL feature columns
N_FEATURES_SUBSET  <- NULL     # analysis: use ALL features
BLOOD_WORKERS      <- 4        # cap PSOCK workers (31GB RAM; unbounded OOMs on 38k-feat VAMS)
here_ <- here::here()

# Protect everything loaded by 0-config.R (run_bioTMLE, Wvars, extract_res, ...)
# plus the control vars -- between-step cleanup only frees the big step data objects.
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
step("prep",       "/src/1 data prep/7-blood-compartment-prep.R")
step("stratified", "/src/2 analysis/12-blood-compartment-intervention-effects.R")
step("combined",   "/src/2 analysis/12-blood-compartment-intervention-effects_combined_arms.R")
step("clean",      "/src/2 analysis/clean_blood_results.R")
cat(sprintf("\n===== ALL_BLOOD_FULL_DONE in %.1f min =====\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
