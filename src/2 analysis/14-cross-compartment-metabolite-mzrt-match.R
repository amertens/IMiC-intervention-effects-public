# =============================================================================
# 14-cross-compartment-metabolite-mzrt-match.R
#
# Cross-compartment concordance of BEP effects on the untargeted METABOLOME.
# Two matching modes:
#   (1) DIRECT ID  -- maternal & infant POSTNATAL VAMS share one `vam_` catalog
#       (same dataset), so they align by feature id, no m/z-RT needed.
#   (2) m/z-RT MATCH -- across matrices with independent Sapient catalogs
#       (milk `rLC` vs plasma `rLC` vs VAMS `vam_`): match on m/z (ppm) + RT.
#
# m/z-RT sources (all available now):
#   maternal plasma / prenatal VAMS  ProcessedDataMISAME3_*.csv (label->MZ,RT)
#   postnatal VAMS (vam_)            metabolite_description_vam_with_global_id.csv (mz,rt_minute,ion mode)
#   milk                             IMiC_alignment.csv (mz_MISAME3 / CHILD_ELICIT_VITAL)  ~68k/97k  [TODO-CONFIRM provenance]
#
# Caveat: cross-PLATFORM RT is only loosely comparable (matrix/column effects);
# m/z is the reliable key, RT a soft filter. Blood is UNADJUSTED for now.
# Outputs: results/cross_compartment_metab_<pair>.csv, figures/cross_compartment/metab_<pair>.png
# =============================================================================

suppressMessages({library(dplyr); library(data.table); library(ggplot2)})
root <- paste0(here::here(), "/"); add <- paste0(root, "data/additional datasets/")
outfig <- paste0(root, "figures/cross_compartment"); dir.create(outfig, showWarnings = FALSE, recursive = TRUE)

# tolerances (CC_PPM=25, method-aware RT, Kim early-elution trim), the m/z-RT
# extractors and match_mzrt() all come from the shared helpers.
source(paste0(root, "src/2 analysis/_blood_helpers.R"))
PPM_TOL <- CC_PPM; RT_TOL_V1V1 <- CC_RT_V1V1; RT_TOL_V1V3 <- CC_RT_V1V3   # local aliases used below

# ---- ATE loaders: feature (UPPERCASED) -> est, sigFDR (representative per feature) ----
ate_milk <- function() {
  readRDS(paste0(root, "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")) %>%
    filter(measure == "ATE", study == "Misame") %>% mutate(feature = toupper(biomarker)) %>%
    group_by(feature) %>% slice_min(pval, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(feature, est, sigFDR)
}
if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- FALSE
.btag <- if (BLOOD_ADJUST) "adjusted_" else ""   # read adjusted_ blood results when set
.osuf <- if (BLOOD_ADJUST) "_adjusted" else ""   # output suffix
ate_blood <- function(dataset, visits) {
  readRDS(paste0(root, "results/blood_compartment_", .btag, "combined_arms_intervention_effects_results_clean.RDS")) %>%
    filter(measure == "ATE", dataset == !!dataset, visit %in% visits) %>% mutate(feature = toupper(biomarker)) %>%
    group_by(feature) %>% slice_min(pval, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(feature, est, sigFDR)
}

# ---- concordance core (shared by m/z-RT and direct-ID pairs) --------------------
report <- function(ov, label) {
  ov <- as.data.table(ov)
  ov[, `:=`(concordant = sign(est_A) == sign(est_B) & est_A != 0 & est_B != 0,
            both_sig = !is.na(sigFDR_A) & !is.na(sigFDR_B) & sigFDR_A == 1 & sigFDR_B == 1)]
  cat(sprintf("[%s] pairs=%d both-FDR-sig=%d concordant(both-sig)=%s\n", label, nrow(ov),
              sum(ov$both_sig), if (sum(ov$both_sig)) sprintf("%.0f%%", 100*mean(ov[both_sig==TRUE]$concordant, na.rm=TRUE)) else "NA"))
  if (nrow(ov)) {
    p <- ggplot(ov, aes(est_A, est_B)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
      geom_point(aes(colour = both_sig), alpha = 0.4, size = 0.8) +
      scale_colour_manual(values = c(`TRUE` = "#E69F00", `FALSE` = "grey75"), guide = "none") +
      labs(x = "ATE (compartment A)", y = "ATE (compartment B)", title = label) + theme_bw()
    ggsave(p, file = paste0(outfig, "/metab_", gsub("[^A-Za-z]", "", label), .osuf, ".png"), width = 7, height = 6.5)
  }
  fwrite(ov, paste0(root, "results/cross_compartment_metab_", gsub("[^A-Za-z]", "", label), .osuf, ".csv"))
  invisible(ov)
}

run_mzrt <- function(mzA, mzB, resA, resB, label, v1v3 = TRUE) {
  pr <- match_mzrt(mzA, mzB, rt_tol = if (v1v3) RT_TOL_V1V3 else RT_TOL_V1V1)
  cat(sprintf("\n=== %s (m/z-RT, %s, %d ppm) ===  matched: %d A-feat -> %d B-feat\n", label,
              if (v1v3) "V1<->V3 loose-RT" else "V1<->V1 tight-RT", PPM_TOL,
              uniqueN(pr$feature_A), uniqueN(pr$feature_B)))
  ov <- as.data.table(pr) |>
    merge(as.data.table(resA)[, .(feature, est_A = est, sigFDR_A = sigFDR)], by.x = "feature_A", by.y = "feature") |>
    merge(as.data.table(resB)[, .(feature, est_B = est, sigFDR_B = sigFDR)], by.x = "feature_B", by.y = "feature")
  report(ov, label)
}
run_id <- function(resA, resB, label) {     # direct id (same catalog)
  cat(sprintf("\n=== %s (direct id) ===\n", label))
  ov <- merge(as.data.table(resA)[, .(feature, est_A = est, sigFDR_A = sigFDR)],
              as.data.table(resB)[, .(feature, est_B = est, sigFDR_B = sigFDR)], by = "feature")
  report(ov, label)
}

# ---- load m/z-RT + ATE once ----------------------------------------------------
milk_mz <- mzrt_milk(misame_only = FALSE); plasma_mz <- mzrt_rlc("ProcessedDataMISAME3_plasma.csv")
prevams_mz <- mzrt_rlc("ProcessedDataMISAME3_VAMS.csv"); infant_mz <- mzrt_vams()
milk_a   <- ate_milk()
plasma_a <- ate_blood("MaternalPlasma", "pn12")
prevam_a <- ate_blood("VamsPrenatal", c("incl","tri3"))
matvam_a <- ate_blood("VamsPostnatalMaternal", "pn56")               # drop tri3 (n=3-9, unstable)
infvam_a <- ate_blood("VamsPostnatalInfant", c("pn12","pn34","pn56"))

# ---- pairs (V1<->V1 = milk/plasma; V1<->V3 = anything <-> postnatal VAMS) -------
run_mzrt(milk_mz, plasma_mz, milk_a, plasma_a, "MilkVsMaternalPlasma", v1v3 = FALSE)  # both LCV1
run_id(matvam_a, infvam_a, "MaternalVAMSvsInfantVAMS")               # same vam_ catalog -> direct
run_mzrt(plasma_mz, infant_mz, plasma_a, infvam_a, "MaternalPlasmaVsInfantVAMS", v1v3 = TRUE)
run_mzrt(milk_mz, infant_mz, milk_a, infvam_a, "MilkVsInfantVAMS", v1v3 = TRUE)       # milk -> infant transfer

cat("\nDone. PPM=25, ionization-mode-matched. RT filtering is OFF by default (match_mzrt\n",
    "use_rt=FALSE, per Trenton's first-pass); the RT_TOL_* args are passed but inactive\n",
    "until use_rt=TRUE. Kim's exact pos/neg RT params still pending.\n")
