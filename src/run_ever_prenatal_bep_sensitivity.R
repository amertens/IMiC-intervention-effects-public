# =============================================================================
# run_ever_prenatal_bep_sensitivity.R
#
# SENSITIVITY: "ever prenatal BEP" contrast for MATERNAL blood compartments.
#   treated  = bep_prenatal == 1  (= BEP/IFA + BEP/BEP)
#   control  = bep_prenatal == 0  (= Control + IFA/BEP)
# Applied to ALL visits (incl. postnatal) -> isolates the PRENATAL axis that the
# period-aware combined coding drops for postnatal maternal samples (where it
# keys on postpartum randomization and buries prenatal-BEP mothers in control).
# Pools to ~n/2 per side, sidestepping the 4-cell small-n instability.
#
# Same machinery as 12-..._combined_arms.R: adjusted Wvars, SL.glm, scale=TRUE,
# bioTMLE by visit, FDR (BH) per visit x dataset. Maternal compartments only.
#
# Out: results/blood_compartment_adjusted_everPrenatalBEP_intervention_effects_results[_clean].RDS
# =============================================================================

rm(list = ls()); source(paste0(here::here(), "/src/0-config.R"))
suppressMessages({library(dplyr); library(data.table)})

BLOOD_WORKERS <- 2                 # avoid the combined-step OOM (38k-feat VAMS)
SCALE_IN_TMLE <- TRUE
Wvars_blood   <- Wvars             # adjusted; Wvars[1]=="arm" is the treatment col
SL_LIB        <- c("SL.glm")

load(file = paste0(here::here(), "/metadata/blood_component.Rdata"))
blood <- readRDS(paste0(here::here(), "/data/blood/merged_blood_datasets.RDS"))
# Postnatal VAMS holds both mother and infant dyad members; keep only mothers below.
vams_postnatal <- blood$vams_postnatal

# MATERNAL compartments only (prenatal randomization is a maternal exposure)
specs <- list(
  MaternalPlasma        = list(df = blood$maternal_plasma,             Y = blood_components$maternal_plasma),
  VamsPrenatal          = list(df = blood$vams_prenatal,               Y = blood_components$vams_prenatal),
  VamsPostnatalMaternal = list(df = dplyr::filter(vams_postnatal, dyad=="mother"), Y = blood_components$vams_postnatal),
  ProteomicsDepleted    = list(df = blood$proteomics_depleted,         Y = blood_components$proteomics_depleted),
  ProteomicsNaive       = list(df = blood$proteomics_naive,            Y = blood_components$proteomics_naive))

run_one <- function(df, Yvars, label) {
  df <- df %>%
    mutate(arm = factor(ifelse(bep_prenatal == 1, "BEP", "Control"),
                        levels = c("Control", "BEP"))) %>%
    filter(!is.na(arm), !is.na(visit))
  if (nrow(df) == 0) { warning(sprintf("[%s] no usable rows", label)); return(NULL) }
  Yvars <- intersect(Yvars, colnames(df))
  cat(sprintf("[%s] n=%d (BEP=%d / Ctrl=%d), %d feats, visits: %s\n",
              label, nrow(df), sum(df$arm=="BEP"), sum(df$arm=="Control"),
              length(Yvars), paste(sort(unique(df$visit)), collapse=",")))
  res <- df %>%
    group_by(visit) %>%
    do(res = try(run_bioTMLE(
      d = ., Wvars = Wvars_blood, bppar.debug = TRUE,
      g_lib = SL_LIB, Q_lib = SL_LIB,
      Yvars = Yvars, scale = SCALE_IN_TMLE,
      bppar.type = BiocParallel::SnowParam(workers = BLOOD_WORKERS))))
  names(res$res) <- paste0(label, "-", res$visit)
  res
}

cat("\n===== ever-prenatal-BEP sensitivity (maternal compartments) =====\n")
results <- lapply(names(specs), function(lab) {
  t <- Sys.time(); r <- run_one(specs[[lab]]$df, specs[[lab]]$Y, lab)
  cat(sprintf("  -> %s done %.1f min\n", lab, as.numeric(difftime(Sys.time(), t, units="mins")))); gc()
  r
})
names(results) <- names(specs)
saveRDS(results, file = paste0(here::here(),
  "/results/blood_compartment_adjusted_everPrenatalBEP_intervention_effects_results.RDS"))

# --- clean: FDR (BH) per visit x dataset (extract_blood_results from helpers) ---
source(paste0(here::here(), "/src/2 analysis/_blood_helpers.R"))
clean <- dplyr::bind_rows(Map(extract_blood_results, results, names(results)))
saveRDS(clean, file = paste0(here::here(),
  "/results/blood_compartment_adjusted_everPrenatalBEP_intervention_effects_results_clean.RDS"))

# --- summary -------------------------------------------------------------------
ate <- clean %>% filter(measure == "ATE")
summ <- ate %>% group_by(dataset, visit) %>%
  summarise(nfeat = n(), nsigFDR = sum(sigFDR==1, na.rm=TRUE),
            pct_up = ifelse(nsigFDR>0, round(100*mean(est[sigFDR==1]>0)), NA_real_),
            .groups="drop") %>% arrange(dataset, visit)
cat("\n===== EVER-PRENATAL-BEP vs NEVER (adjusted), maternal compartments =====\n")
print(as.data.frame(summ))
cat("\nSaved: results/blood_compartment_adjusted_everPrenatalBEP_intervention_effects_results_clean.RDS\n")
