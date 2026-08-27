# =============================================================================
# 19-cross-compartment-arrow-contrast-summary.R   [SUPPORTING — demoted 2026-06-25]
#
# NOTE: superseded as the primary by 19b (FDR-first annotated lists). This bulk
# %-concordance table is kept as supporting context. Matching is now ppm-only
# (25 ppm + ionization mode, no RT) via the shared helper.
#
# Trenton's "hard-code this table" (Kim x Trenton meeting): for each FORWARD
# temporality arrow (maternal blood -> milk -> infant blood) x each group
# contrast, count significant features (p<0.05, NOT FDR) per side, how many
# match within 25 ppm (+ ionization mode + method-aware RT), and the % of
# matched pairs that are direction-concordant.
#
# Contrasts (stratified arms vs Control): BEP/BEP, IFA/BEP (postnatal-only),
# BEP/IFA (prenatal-only). Uses the milk + blood STRATIFIED clean results.
# Matching reuses the 25-ppm / mode / RT-trim rules from script 14.
#
# Out: results/cross_compartment_arrow_contrast_summary.csv
#
# NOTE (2026-06-25): this bulk p<0.05 table is now the DEMOTED supporting check.
# The FDR-FIRST per-contrast annotated lists (the primary deliverable per Trenton)
# are produced by 19b-cross-compartment-fdr-first-lists.R. Matching is ppm-only
# (use_rt=FALSE in _blood_helpers) per Trenton's "wouldn't worry about RT" steer.
# =============================================================================

suppressMessages({library(data.table)})
root <- paste0(here::here(), "/"); add <- paste0(root, "data/additional datasets/")

# tolerances, m/z-RT extractors and the matcher all live in the shared helpers
source(paste0(root, "src/2 analysis/_blood_helpers.R"))
RT_TOL_V1V1 <- CC_RT_V1V1; RT_TOL_V1V3 <- CC_RT_V1V3   # local aliases used below

# --- stratified results (carry the 4-level `contrast`) --------------------------
milk  <- as.data.table(readRDS(paste0(root, "results/adjusted_intervention_effects_res_untargeted_metabolomics_clean.RDS")))
blood <- as.data.table(readRDS(paste0(root, "results/blood_compartment_adjusted_intervention_effects_results_clean.RDS")))
milk  <- milk[measure == "ATE" & study == "Misame"]
blood <- blood[measure == "ATE"]

# Feature-level table for one milk contrast x visit. When a feature appears more
# than once, keep its most significant (smallest-pval) row as the representative.
get_milk <- function(cc, vv) {
  d <- milk[contrast == cc & visit == vv]
  d[, feature := toupper(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]
  d[, .(feature, est, pval, sigFDR)]
}
# Same, for one blood dataset x visit(s) x contrast.
get_blood <- function(dd, vv, cc) {
  d <- blood[dataset == dd & visit %in% vv & contrast == cc]
  d[, feature := toupper(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]
  d[, .(feature, est, pval, sigFDR)]
}

# m/z-RT tables
mk_mz <- mzrt_milk(misame_only = FALSE); pl_mz <- mzrt_rlc("ProcessedDataMISAME3_plasma.csv"); vm_mz <- mzrt_vams()
pv_mz <- mzrt_rlc("ProcessedDataMISAME3_VAMS.csv")   # maternal PRENATAL VAMS (V1)

# --- arrows: forward temporality (maternal blood -> milk -> infant blood) --------
# temp = contemporaneous (same calendar window) | forward-lag (upstream earlier).
# Age key: pn12=1-2mo, pn34=3-4mo, pn56=5-6mo postpartum; milk "1-2 mo."/"3-4 mo.".
# Data limit: maternal blood postpartum exists only at pn12 (1-2mo) for plasma.
arrows <- list(
  # ---- CONTEMPORANEOUS (same window) ----
  list(name="Plasma pn12 -> Milk 1-2mo",     temp="contemporaneous", A=list("blood","MaternalPlasma","pn12",pl_mz),        B=list("milk",NA,"1-2 mo.",mk_mz),                   match="v1v1"),
  list(name="Milk 1-2mo -> Infant pn12",     temp="contemporaneous", A=list("milk",NA,"1-2 mo.",mk_mz),                   B=list("blood","VamsPostnatalInfant","pn12",vm_mz),  match="v1v3"),
  list(name="Plasma pn12 -> Infant pn12",    temp="contemporaneous", A=list("blood","MaternalPlasma","pn12",pl_mz),        B=list("blood","VamsPostnatalInfant","pn12",vm_mz),  match="v1v3"),
  list(name="Milk 3-4mo -> Infant pn34",     temp="contemporaneous", A=list("milk",NA,"3-4 mo.",mk_mz),                   B=list("blood","VamsPostnatalInfant","pn34",vm_mz),  match="v1v3"),
  list(name="MatVAMS pn56 <-> Infant pn56",  temp="contemporaneous", A=list("blood","VamsPostnatalMaternal","pn56",vm_mz), B=list("blood","VamsPostnatalInfant","pn56",vm_mz),  match="id"),
  # ---- FORWARD-LAG (upstream earlier than downstream) ----
  list(name="Plasma pn12 -> Milk 3-4mo",     temp="forward-lag",     A=list("blood","MaternalPlasma","pn12",pl_mz),        B=list("milk",NA,"3-4 mo.",mk_mz),                   match="v1v1"),
  list(name="Milk 1-2mo -> Infant pn34",     temp="forward-lag",     A=list("milk",NA,"1-2 mo.",mk_mz),                   B=list("blood","VamsPostnatalInfant","pn34",vm_mz),  match="v1v3"),
  list(name="Milk 3-4mo -> Infant pn56",     temp="forward-lag",     A=list("milk",NA,"3-4 mo.",mk_mz),                   B=list("blood","VamsPostnatalInfant","pn56",vm_mz),  match="v1v3"),
  list(name="Plasma pn12 -> Infant pn34",    temp="forward-lag",     A=list("blood","MaternalPlasma","pn12",pl_mz),        B=list("blood","VamsPostnatalInfant","pn34",vm_mz),  match="v1v3"),
  list(name="Plasma pn12 -> Infant pn56",    temp="forward-lag",     A=list("blood","MaternalPlasma","pn12",pl_mz),        B=list("blood","VamsPostnatalInfant","pn56",vm_mz),  match="v1v3"),
  # ---- PRENATAL-LAG (maternal prenatal blood -> infant birth; matches the mock Tri2->Birth example).
  #      Maternal prenatal blood is V1, infant birth VAMS is V3 -> cross-platform mass-match (no shared catalogue).
  #      At birth only prenatal-BEP contrasts (BEP/BEP, BEP/IFA) are interpretable; IFA/BEP is a negative control. ----
  list(name="MatVAMS incl -> Infant acco",   temp="prenatal-lag",    A=list("blood","VamsPrenatal","incl",pv_mz),          B=list("blood","VamsPostnatalInfant","acco",vm_mz), match="v1v3"),
  list(name="MatVAMS tri3 -> Infant acco",   temp="prenatal-lag",    A=list("blood","VamsPrenatal","tri3",pv_mz),          B=list("blood","VamsPostnatalInfant","acco",vm_mz), match="v1v3"),
  list(name="Plasma tri3 -> Infant acco",    temp="prenatal-lag",    A=list("blood","MaternalPlasma","tri3",pl_mz),        B=list("blood","VamsPostnatalInfant","acco",vm_mz), match="v1v3")
)
CONTRASTS <- c("BEP/BEP", "IFA/BEP", "BEP/IFA")

side_feats <- function(spec, cc) {
  if (spec[[1]] == "milk") get_milk(cc, spec[[3]]) else get_blood(spec[[2]], spec[[3]], cc)
}

# match two feature-sets, return (npair, nAm, pconc) -- mtype id/v1v1/v1v3
match_pairs <- function(selA, selB, mtype, mzA_tab, mzB_tab) {
  if (nrow(selA) == 0 || nrow(selB) == 0) return(list(npair=0L, nAm=0L, pconc=NA_real_))
  if (mtype == "id") {
    m <- merge(selA, selB, by = "feature", suffixes = c("_A","_B"))
    if (!nrow(m)) return(list(npair=0L, nAm=0L, pconc=NA_real_))
    return(list(npair=nrow(m), nAm=nrow(m), pconc=round(100*mean(sign(m$est_A) == sign(m$est_B)))))
  }
  rt <- if (mtype == "v1v1") RT_TOL_V1V1 else RT_TOL_V1V3
  pr <- match_mzrt(mzA_tab[feature %in% selA$feature], mzB_tab[feature %in% selB$feature], rt)
  if (!nrow(pr)) return(list(npair=0L, nAm=0L, pconc=NA_real_))
  pr <- merge(pr, selA[, .(feature, est_A = est)], by.x = "feature_A", by.y = "feature")
  pr <- merge(pr, selB[, .(feature, est_B = est)], by.x = "feature_B", by.y = "feature")
  list(npair=nrow(pr), nAm=uniqueN(pr$feature_A), pconc=round(100*mean(sign(pr$est_A) == sign(pr$est_B))))
}

rows <- list()
for (ar in arrows) for (cc in CONTRASTS) {
  fA <- side_feats(ar$A, cc); fB <- side_feats(ar$B, cc)
  p05 <- match_pairs(fA[pval < 0.05], fB[pval < 0.05], ar$match, ar$A[[4]], ar$B[[4]])     # noisy bulk
  fdr <- match_pairs(fA[sigFDR == 1], fB[sigFDR == 1], ar$match, ar$A[[4]], ar$B[[4]])     # real-signal
  nA05 <- nrow(fA[pval < 0.05])
  rows[[length(rows)+1]] <- data.table(temporality=ar$temp, arrow=ar$name, contrast=cc, match=ar$match,
    n_sig_A=nA05, n_sig_B=nrow(fB[pval < 0.05]),
    n_matched_pairs=p05$npair, pct_A_matched=ifelse(nA05>0, round(100*p05$nAm/nA05), NA_real_),
    pct_concordant_p05=p05$pconc,
    n_FDR_A=nrow(fA[sigFDR == 1]), n_FDR_B=nrow(fB[sigFDR == 1]),
    n_matched_FDR=fdr$npair, pct_concordant_FDR=fdr$pconc)
}
tab <- rbindlist(rows)
fwrite(tab, paste0(root, "results/cross_compartment_arrow_contrast_summary.csv"))
cat("=== Cross-compartment arrow x contrast summary (25 ppm) ===\n")
print(tab[], nrow = 100)
cat("\npct_concordant_p05 = direction concordance over the p<0.05 matched bulk (noise-limited for cross-platform);\n")
cat("pct_concordant_FDR = same over the FDR-significant subset (the real-signal version).\n")
cat("Saved: results/cross_compartment_arrow_contrast_summary.csv\n")
