# =============================================================================
# 39-four-compartment-bep-feature-screen.R
#
# Consolidated "BEP product -> maternal blood -> milk -> infant blood" feature screen.
#
# REVIEWER-2 mechanistic ask (2026-06-30 IMiC meeting, Trenton): find a metabolite
# that is ABUNDANT IN THE BEP PRODUCT *and* elevated in INFANT BLOOD, with milk and
# maternal blood as the intervening compartments. In a randomized trial that chain is
# essentially irrefutable -- "no question it got there via mom's breast milk".
#
# The pairwise legs already exist (14 = cross-compartment m/z-RT concordance;
# 21 = supplement-tracer enrichment; 22 = supplement-composition cross-ref). This
# script ties them into ONE ranked scorecard, anchored on each infant-blood (V3 vam_)
# feature that BEP RAISED, attaching for each:
#   up_milk     BEP raised the m/z-matched milk feature           (V1<->V3 isobaric)
#   up_plasma   BEP raised the m/z-matched maternal-plasma feature (V1<->V3 isobaric)
#   up_matVAMS  BEP raised the SAME vam_ feature in the mother      (V3 identity, cleanest)
#   in_product  an m/z-matched feature is present/abundant in the BEP product (V1<->V3)
# full_chain  = up in infant AND (plasma OR matVAMS) AND milk AND present in product.
#
# Matching is ppm-only (25 ppm + ionization mode), per Trenton "don't worry about RT
# yet". Cross-platform (V1<->V3) legs are ISOBARIC-LIMITED, so the identity
# maternal-VAMS leg and any FDR-sig hits are the high-confidence evidence; the rest is
# a candidate list to annotate. Direction-only (Kim intensity constraint). MISAME-III,
# BEP-vs-control combined-arms (adjusted by default).
#
# EXPECTED RESULT (consistent with script 22): most BEP-up features are NOT abundant
# product components -> a short/empty full-chain list is itself the honest answer
# (endogenous metabolic response, e.g. the octenoylcarnitine pair). The selenium /
# selenoprotein chain remains the cleaner ready-made example (proteome, script 13).
#
# Out: results/four_compartment_bep_feature_screen.csv  (ranked scorecard, all infant-up)
#      results/four_compartment_bep_shortlist.csv        (features passing the full chain)
# =============================================================================

suppressMessages({library(data.table)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))

if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- TRUE     # primary analysis = adjusted
.btag <- if (BLOOD_ADJUST) "adjusted_" else ""

# -- thresholds ----------------------------------------------------------------
PRODUCT_PRESENT <- 0     # supp_mean >= 0  -> at/above maternal-VAMS background ("present")
PRODUCT_ABUND   <- 1     # supp_mean >= 1  -> >=1 SD above background ("abundant in product")

# -- ATE loaders (one representative row per feature: smallest p) ---------------
# A feature can appear at several visits; keep only its single most-significant row
# so each feature contributes one direction/effect to the scorecard.
rep_feat <- function(d) {
  d <- as.data.table(d)
  d[, feature := toupper(biomarker)]
  d <- d[order(pval)]                 # smallest p first
  d[!duplicated(feature)]             # keep the first (most-significant) row per feature
}

# Blood ATEs: BEP-vs-control contrast only, one file spanning all blood compartments.
blood <- as.data.table(readRDS(paste0(root, "results/blood_compartment_", .btag,
            "combined_arms_intervention_effects_results_clean.RDS")))[measure == "ATE" & contrast == "BEP"]
inf  <- rep_feat(blood[dataset == "VamsPostnatalInfant"   & visit %in% c("pn12", "pn34", "pn56")])  # infant blood (anchor)
matv <- rep_feat(blood[dataset == "VamsPostnatalMaternal" & visit == "pn56"])                        # maternal VAMS (identity leg)
plas <- rep_feat(blood[dataset == "MaternalPlasma"        & visit == "pn12"])                        # maternal plasma
milk <- rep_feat(as.data.table(readRDS(paste0(root,
            "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")))[
            measure == "ATE" & study == "Misame" & contrast == "BEP" & visit %in% c("1-2 mo.", "3-4 mo.")])

# -- m/z-RT catalogues + supplement composition --------------------------------
inf_mz  <- mzrt_vams()
plas_mz <- mzrt_rlc("ProcessedDataMISAME3_plasma.csv")
milk_mz <- mzrt_milk(misame_only = FALSE)
supp    <- fread(paste0(root, "results/bep_supplement_profile_V1.csv")); supp[, feature := toupper(feature)]
supp_mz <- merge(mzrt_rlc("ProcessedDataMISAME3_VAMS.csv"), supp[, .(feature, supp_mean)], by = "feature")

# -- anchor: infant features BEP RAISED (est > 0); carry the FDR flag ----------
anchor    <- inf[est > 0]
anchor_mz <- inf_mz[feature %in% anchor$feature]
cat(sprintf("Infant BEP-up features: %d (FDR-sig %d) | with m/z-RT: %d\n",
            nrow(anchor), sum(anchor$sigFDR == 1), nrow(anchor_mz)))

# -- attach a V1 compartment to the infant anchor via best 25-ppm match --------
# For each infant anchor feature, find the nearest-ppm upstream-compartment feature
# (best_match returns feature_A = infant, feature_B = upstream, 1:1) and carry that
# upstream feature's effect + FDR flag back onto the infant feature. `pfx` namespaces
# the columns (e.g. milk_est / milk_sigFDR) so several compartments can be attached.
attach_v1 <- function(tgt_res, tgt_mz, pfx) {
  matched_pairs <- best_match(anchor_mz, tgt_mz[feature %in% tgt_res$feature])
  if (!nrow(matched_pairs)) {
    out <- data.table(feature = character(), est = numeric(), cil = numeric(),
                      ciu = numeric(), sigFDR = numeric())
  } else {
    upstream_effects <- tgt_res[, .(feature_B = feature, est, cil, ciu, sigFDR)]
    out <- merge(matched_pairs, upstream_effects, by = "feature_B")[
             , .(feature = feature_A, est, cil, ciu, sigFDR)]  # relabel back to the infant feature id
  }
  setnames(out, c("est", "cil", "ciu", "sigFDR"), paste0(pfx, c("_est", "_cil", "_ciu", "_sigFDR")))
  out
}
# supplement abundance attached the same way (V1<->V3 isobaric)
attach_supp <- function() {
  matched_pairs <- best_match(anchor_mz, supp_mz)
  if (!nrow(matched_pairs)) return(data.table(feature = character(), supp_mean = numeric()))
  merge(matched_pairs, supp_mz[, .(feature_B = feature, supp_mean)], by = "feature_B")[
    , .(feature = feature_A, supp_mean)]
}

# `sc` = the scorecard: one row per infant anchor feature, gaining an upstream-
# compartment column block with each merge (all.x keeps every anchor feature).
sc <- anchor[, .(feature, infant_visit = visit, infant_est = est, infant_cil = cil,
                 infant_ciu = ciu, infant_sigFDR = sigFDR)]
sc <- merge(sc, attach_v1(milk, milk_mz, "milk"),   by = "feature", all.x = TRUE)
sc <- merge(sc, attach_v1(plas, plas_mz, "plasma"), by = "feature", all.x = TRUE)
sc <- merge(sc, matv[, .(feature, matVAMS_est = est, matVAMS_cil = cil, matVAMS_ciu = ciu,
                         matVAMS_sigFDR = sigFDR)], by = "feature", all.x = TRUE)  # V3 identity (no ppm match needed)
sc <- merge(sc, attach_supp(), by = "feature", all.x = TRUE)

# -- scorecard flags -----------------------------------------------------------
# Each leg is "raised by BEP" only when it matched (est not NA) AND the effect is positive.
sc[, up_milk    := !is.na(milk_est)    & milk_est    > 0]
sc[, up_plasma  := !is.na(plasma_est)  & plasma_est  > 0]
sc[, up_matVAMS := !is.na(matVAMS_est) & matVAMS_est > 0]
sc[, up_maternal := up_plasma | up_matVAMS]
sc[, in_product           := !is.na(supp_mean) & supp_mean >= PRODUCT_PRESENT]
sc[, abundant_in_product  := !is.na(supp_mean) & supp_mean >= PRODUCT_ABUND]
sc[, n_criteria := up_maternal + up_milk + in_product]                  # infant-up is the anchor (always TRUE)
sc[, full_chain          := up_maternal & up_milk & in_product]
sc[, full_chain_abundant := up_maternal & up_milk & abundant_in_product]

# -- annotate infant feature (m/z, RT, compound class/name where available) ----
desc <- fread(paste0(root, "data/additional datasets/metabolite_description_vam_with_global_id.csv"))
desc[, feature := toupper(feature_label)]
namecols <- intersect(c("Compound Name", "compound_name", "name", "Metabolite",
                        "Compound Class", "global_id"), names(desc))
sc <- merge(sc, desc[, c("feature", "mz", "rt_minute", namecols), with = FALSE], by = "feature", all.x = TRUE)

setorder(sc, -full_chain_abundant, -full_chain, -n_criteria, -infant_sigFDR, -supp_mean, na.last = TRUE)
fwrite(sc, paste0(root, "results/four_compartment_bep_feature_screen.csv"))
short <- sc[full_chain == TRUE]
fwrite(short, paste0(root, "results/four_compartment_bep_shortlist.csv"))

# -- summary -------------------------------------------------------------------
cat("\n=== Four-compartment BEP feature screen (anchored on infant BEP-up features) ===\n")
cat(sprintf("  up in milk (matched)            : %d\n", sum(sc$up_milk)))
cat(sprintf("  up in maternal (plasma|VAMS)    : %d  (identity matVAMS: %d)\n", sum(sc$up_maternal), sum(sc$up_matVAMS)))
cat(sprintf("  present in product (supp>=%g)    : %d\n", PRODUCT_PRESENT, sum(sc$in_product)))
cat(sprintf("  abundant in product (supp>=%g)   : %d\n", PRODUCT_ABUND, sum(sc$abundant_in_product)))
cat(sprintf("  FULL CHAIN (present in product)  : %d\n", sum(sc$full_chain)))
cat(sprintf("  FULL CHAIN (ABUNDANT in product) : %d   <- the irrefutable-linkage candidates\n", sum(sc$full_chain_abundant)))
cat(sprintf("  ...of which ROBUST (infant FDR-sig): chain=%d  abundant-chain=%d   <- decision-relevant\n",
            sum(sc$full_chain & sc$infant_sigFDR == 1), sum(sc$full_chain_abundant & sc$infant_sigFDR == 1)))
cat("\nTop full-chain candidates:\n")
print(head(sc[full_chain == TRUE, .(feature, mz, infant_sigFDR, up_matVAMS, supp_mean,
                                    abundant_in_product)], 15))
cat("\nOcteno(y)l-carnitine reference (infant vam_1005523 / vam_1005524, m/z ~286.20):\n")
print(sc[feature %in% c("VAM_1005523", "VAM_1005524"),
         .(feature, infant_est, infant_sigFDR, up_milk, up_matVAMS, up_plasma, supp_mean, in_product)])
cat("\nSaved: results/four_compartment_bep_feature_screen.csv (+ _shortlist.csv)\n")
