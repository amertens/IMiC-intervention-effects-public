# =============================================================================
# 20-cross-compartment-threshold-free-panel.R   [SUPPORTING, ppm-only 2026-06-25]
#
# Supporting (not the primary; see 19b). Matching is ppm-only (25 ppm + mode, no RT).
#
# Threshold-free concordance (RRHO + GSEA + weighted r + anchored) across the
# SAME forward-temporality arrows x contrasts as the crude table (script 19),
# so Trenton sees noise-robust numbers next to the crude % concordant.
#
# Per arrow x contrast, over matched feature PAIRS (id = shared vam_ catalogue;
# cross-platform = best 1:1 m/z-RT match at 25 ppm + mode + RT), using the FULL
# ranked effect distributions (no FDR threshold needed):
#   - weighted Pearson r (+ Deming slope)
#   - RRHO directional: concordant vs discordant peak (-log10 hypergeom p)
#   - GSEA: A's nominal-up signature enriched in B's signed ranking (perm p/NES)
#   - anchored sign test among A's p<0.05 features
#
# Out: results/cross_compartment_threshold_free_panel.csv
#
# NOTE (2026-06-25): supporting (demoted) full-distribution panel; the FDR-FIRST
# per-contrast annotated lists (primary deliverable) are in 19b-...-fdr-first-lists.R.
# Matching is ppm-only (use_rt=FALSE in _blood_helpers) per Trenton's RT steer.
# =============================================================================

suppressMessages({library(data.table)})
root <- paste0(here::here(), "/"); add <- paste0(root, "data/additional datasets/")
# tolerances, m/z-RT extractors, best_match, and the threshold-free methods
# (signed_stat, gsea_es, perm_gsea_p, rrho_directional, weighted_cor, deming_slope)
# all come from the shared helpers.
source(paste0(root, "src/2 analysis/_blood_helpers.R"))
RT_V1V1 <- CC_RT_V1V1; RT_V1V3 <- CC_RT_V1V3   # local aliases used in the arrow loop
set.seed(20240624)                              # reproducible GSEA permutation p-values

# threshold-free method panel for one paired table (est/pval/se per side)
methods <- function(dt){
  x <- dt$est_A; y <- dt$est_B; n <- nrow(dt)
  # inverse-variance weights combining both sides' SEs; drop any non-finite weight
  w <- 1/(dt$se_A^2 + dt$se_B^2); w[!is.finite(w)] <- NA

  weighted_r <- weighted_cor(x, y, w)
  deming     <- deming_slope(x, y)
  rrho       <- rrho_directional(signed_stat(dt$pval_A, dt$est_A), signed_stat(dt$pval_B, dt$est_B))
  # GSEA: is A's nominal-up signature (p<0.05 & est>0) enriched in B's signed ranking?
  gsea       <- perm_gsea_p(signed_stat(dt$pval_B, dt$est_B), dt$pval_A < 0.05 & dt$est_A > 0, P = 300)

  # Anchored sign test: among A's p<0.05 features, how often does B agree in direction?
  anchor <- dt$pval_A < 0.05
  if (sum(anchor) > 2) {
    n_concordant <- sum(sign(x[anchor]) == sign(y[anchor]))
    anchored_pct <- round(100 * mean(sign(x[anchor]) == sign(y[anchor])))
    anchored_p   <- binom.test(n_concordant, sum(anchor), 0.5)$p.value
  } else {
    anchored_pct <- NA; anchored_p <- NA
  }

  data.table(n_pairs = n, w_r = round(weighted_r, 3), deming = round(deming, 2),
             rrho_conc = round(rrho[["concordant"]], 1), rrho_disc = round(rrho[["discordant"]], 1),
             gsea_up_NES = round(gsea$nes, 2), gsea_up_p = signif(gsea$p, 2),
             anchored_pct = anchored_pct, anchored_p = signif(anchored_p, 2))
}

# --- results + extractors ------------------------------------------------------
milk  <- as.data.table(readRDS(paste0(root,"results/adjusted_intervention_effects_res_untargeted_metabolomics_clean.RDS")))[measure=="ATE" & study=="Misame"]
blood <- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_intervention_effects_results_clean.RDS")))[measure=="ATE"]
# standard error backed out of a 95% CI half-width
se_of <- function(cil,ciu) (ciu-cil)/(2*qnorm(0.975))
# Feature-level tables (representative = smallest-pval row per feature), carrying an SE.
get_milk  <- function(cc,vv){
  d <- milk[contrast==cc & visit==vv]
  d[, feature := toupper(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]
  d[, .(feature, est, pval, se=se_of(cil,ciu))]
}
get_blood <- function(dd,vv,cc){
  d <- blood[dataset==dd & visit %in% vv & contrast==cc]
  d[, feature := toupper(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]
  d[, .(feature, est, pval, se=se_of(cil,ciu))]
}
mk_mz<-mzrt_milk(misame_only=FALSE); pl_mz<-mzrt_rlc("ProcessedDataMISAME3_plasma.csv"); vm_mz<-mzrt_vams()

arrows <- list(
  # ---- CONTEMPORANEOUS (same window) ----
  list(name="Plasma pn12 -> Milk 1-2mo",     temp="contemporaneous", A=list("blood","MaternalPlasma","pn12",pl_mz),        B=list("milk",NA,"1-2 mo.",mk_mz),                  match="v1v1"),
  list(name="Milk 1-2mo -> Infant pn12",     temp="contemporaneous", A=list("milk",NA,"1-2 mo.",mk_mz),                   B=list("blood","VamsPostnatalInfant","pn12",vm_mz), match="v1v3"),
  list(name="Plasma pn12 -> Infant pn12",    temp="contemporaneous", A=list("blood","MaternalPlasma","pn12",pl_mz),        B=list("blood","VamsPostnatalInfant","pn12",vm_mz), match="v1v3"),
  list(name="Milk 3-4mo -> Infant pn34",     temp="contemporaneous", A=list("milk",NA,"3-4 mo.",mk_mz),                   B=list("blood","VamsPostnatalInfant","pn34",vm_mz), match="v1v3"),
  list(name="MatVAMS pn56 <-> Infant pn56",  temp="contemporaneous", A=list("blood","VamsPostnatalMaternal","pn56",vm_mz), B=list("blood","VamsPostnatalInfant","pn56",vm_mz), match="id"),
  # ---- FORWARD-LAG (upstream earlier than downstream) ----
  list(name="Plasma pn12 -> Milk 3-4mo",     temp="forward-lag",     A=list("blood","MaternalPlasma","pn12",pl_mz),        B=list("milk",NA,"3-4 mo.",mk_mz),                  match="v1v1"),
  list(name="Milk 1-2mo -> Infant pn34",     temp="forward-lag",     A=list("milk",NA,"1-2 mo.",mk_mz),                   B=list("blood","VamsPostnatalInfant","pn34",vm_mz), match="v1v3"),
  list(name="Milk 3-4mo -> Infant pn56",     temp="forward-lag",     A=list("milk",NA,"3-4 mo.",mk_mz),                   B=list("blood","VamsPostnatalInfant","pn56",vm_mz), match="v1v3"),
  list(name="Plasma pn12 -> Infant pn34",    temp="forward-lag",     A=list("blood","MaternalPlasma","pn12",pl_mz),        B=list("blood","VamsPostnatalInfant","pn34",vm_mz), match="v1v3"),
  list(name="Plasma pn12 -> Infant pn56",    temp="forward-lag",     A=list("blood","MaternalPlasma","pn12",pl_mz),        B=list("blood","VamsPostnatalInfant","pn56",vm_mz), match="v1v3")
)
CONTRASTS <- c("BEP/BEP","IFA/BEP","BEP/IFA")
gf <- function(spec,cc) if(spec[[1]]=="milk") get_milk(cc,spec[[3]]) else get_blood(spec[[2]],spec[[3]],cc)

out <- list()
for (ar in arrows) {
  # Build the matched-pair skeleton ONCE per arrow (matching is contrast-independent):
  # id arrows share the vam_ catalogue; cross-platform arrows use a 1:1 m/z-RT match.
  fa0 <- gf(ar$A, CONTRASTS[1]); fb0 <- gf(ar$B, CONTRASTS[1])
  if (ar$match == "id") {
    skel <- data.table(feature_A = intersect(fa0$feature, fb0$feature))[, feature_B := feature_A]
  } else {
    rt   <- if (ar$match == "v1v1") RT_V1V1 else RT_V1V3
    skel <- best_match(ar$A[[4]][feature %in% fa0$feature], ar$B[[4]][feature %in% fb0$feature], rt)
  }
  for (cc in CONTRASTS) {
    # attach each side's per-contrast effect / SE onto the fixed skeleton
    fA <- gf(ar$A, cc); fB <- gf(ar$B, cc)
    dt <- merge(skel, fA[, .(feature_A = feature, est_A = est, pval_A = pval, se_A = se)], by = "feature_A")
    dt <- merge(dt,   fB[, .(feature_B = feature, est_B = est, pval_B = pval, se_B = se)], by = "feature_B")
    dt <- dt[is.finite(est_A) & is.finite(est_B) & is.finite(pval_A) & is.finite(pval_B)]
    if (nrow(dt) < 20) {   # too few pairs for stable threshold-free stats -> emit an all-NA row
      out[[length(out)+1]] <- data.table(temporality=ar$temp, arrow=ar$name, contrast=cc, n_pairs=nrow(dt),
                                          w_r=NA, deming=NA, rrho_conc=NA, rrho_disc=NA,
                                          gsea_up_NES=NA, gsea_up_p=NA, anchored_pct=NA, anchored_p=NA)
      next
    }
    out[[length(out)+1]] <- cbind(temporality=ar$temp, arrow=ar$name, contrast=cc, methods(dt))
  }
  cat("done:", ar$name, "\n")
}
panel <- rbindlist(out)
fwrite(panel, paste0(root,"results/cross_compartment_threshold_free_panel.csv"))
cat("\n=== Threshold-free panel (full ranked distributions over matched pairs) ===\n")
print(panel[], nrow=100)
cat("\nrrho_conc/disc = -log10 hypergeom p (concordant vs discordant corner); gsea_up = A-up signature in B ranking;\n")
cat("anchored = sign concordance among A p<0.05. Saved: results/cross_compartment_threshold_free_panel.csv\n")
