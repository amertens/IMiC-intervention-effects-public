# =============================================================================
# 22-supplement-composition-crossref.R
#
# PROVENANCE / TRACER question (Andrew, for Trenton): are the features the
# intervention INCREASED in mother & infant (postnatal VAMS, V3) the ones that are
# ABUNDANT in the BEP supplement product itself (results/bep_supplement_profile_V1.csv,
# V1)? If yes -> direct dietary transfer; if not -> endogenous metabolic response.
#
# Platform: supplement = V1 (RLC_); postnatal-VAMS up-features = V3 (vam_) -> the
# cross-reference is CROSS-PLATFORM accurate-mass (25 ppm + ionization mode),
# ISOBARIC-LIMITED. supp_mean = z-scored intensity in product vs maternal VAMS
# (high = abundant in product). MISAME-III only; direction-only (Kim constraint).
#
# Out: results/bep_supplement_composition_xref.csv  (per up-feature: best supp match)
# =============================================================================

suppressMessages({library(data.table)})
root <- paste0(here::here(), "/"); source(paste0(root, "src/2 analysis/_blood_helpers.R"))

# -- supplement composition (V1) + m/z-RT from the prenatal-VAMS V1 catalogue ----
supp <- fread(paste0(root, "results/bep_supplement_profile_V1.csv")); supp[, feature := toupper(feature)]
N <- nrow(supp)
supp <- merge(supp, mzrt_rlc("ProcessedDataMISAME3_VAMS.csv"), by = "feature", all.x = TRUE)
supp[, pct_top := 100*(N - supp_rank + 1)/N]                  # 100 = most abundant in product
suppM <- supp[!is.na(mz) & !is.na(mode)]

# -- intervention-UP V3 features (arm-stratified, FDR-sig & est>0) ---------------
up_ate <- as.data.table(readRDS(paste0(root, "results/blood_compartment_adjusted_intervention_effects_results_clean.RDS")))[measure=="ATE" & sigFDR==1 & est>0]
vm <- mzrt_vams()
up_set <- function(ds, vis, lbl) {
  f <- toupper(unique(up_ate[dataset==ds & visit %in% vis]$biomarker))
  cbind(set = lbl, merge(data.table(feature=f), vm, by="feature", all.x=TRUE))
}
ups <- rbind(up_set("VamsPostnatalMaternal","pn56","Maternal pn56 UP"),
             up_set("VamsPostnatalInfant",c("pn12","pn34","pn56"),"Infant UP"))

# -- best cross-platform supplement match (25 ppm + mode) per up-feature ---------
best_supp <- function(mz0, mode0) {
  no_match <- data.table(n_supp_match=0L, best_supp_mean=NA_real_, best_pct_top=NA_real_)
  if (is.na(mz0) || is.na(mode0)) return(no_match)
  # supplement features in the same ionization mode within 25 ppm of this m/z
  cand <- suppM[mode==mode0 & abs(mz-mz0)/mz0*1e6 <= CC_PPM]
  if (!nrow(cand)) return(no_match)
  cand[, .(n_supp_match=.N, best_supp_mean=max(supp_mean), best_pct_top=max(pct_top))]
}
ups <- cbind(ups, rbindlist(Map(best_supp, ups$mz, ups$mode)))
ups[, in_supp_high := best_supp_mean >= 1]
fwrite(ups[order(set, -best_supp_mean)], paste0(root, "results/bep_supplement_composition_xref.csv"))

# -- verdict summary ------------------------------------------------------------
cat("=== Supplement composition cross-reference (V3 up-features -> V1 supplement, 25 ppm) ===\n")
ups[, .(n=.N, high_in_supp=sum(best_supp_mean>=1,na.rm=TRUE),
        no_match=sum(n_supp_match==0|is.na(best_supp_mean)),
        median_best_supp_mean=round(median(best_supp_mean,na.rm=TRUE),2)), by=set][] |> print()
cat("\nCarnitine m/z 286.202 (pos) supplement matches within 25 ppm:",
    nrow(suppM[mode=="positive" & abs(mz-286.202)/286.202*1e6<=CC_PPM]), "\n")
cat("Saved: results/bep_supplement_composition_xref.csv\n")
