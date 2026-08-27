# =============================================================================
# 28-matching-sensitivity.R
#
# Does better feature linkage firm up the cross-compartment agreement? Two cheap
# precision upgrades over the mass-only (25 ppm) first pass, on the combined-arms
# (any-postnatal-BEP) effects:
#   (1) SAME-PLATFORM pair (milk <-> maternal plasma, both V1): add the RT window
#       back (valid because same chromatography) and compare to mass-only.
#   (2) CROSS-PLATFORM pairs (milk/plasma <-> infant, V1<->V3, RT not comparable):
#       restrict to UNIQUELY-matchable masses (exactly one candidate within 25 ppm
#       on BOTH sides) and compare to the full mass-only match.
# Metric: % of matched features agreeing in BEP direction (all, and among the
# features significant in A), plus weighted effect correlation.
#
# Out: results/cross_compartment_matching_sensitivity.csv
# =============================================================================
suppressMessages({library(data.table)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))

# Feature IDs are matched case-insensitively across catalogues; upper-case is the key.
feature_key <- function(x) toupper(as.character(x))
# Standard error implied by a symmetric 95% CI (est +/- 1.96*SE).
se_of <- function(cil, ciu) (ciu - cil)/(2*qnorm(0.975))

milkC  <- as.data.table(readRDS(paste0(root,"results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")))
bloodC <- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))

# Pull one milk visit's ATE effects as one row per feature, keeping the most
# significant duplicate. Returns feature key, effect, p-value and implied SE.
get_milk <- function(vv){
  d <- milkC[measure=="ATE" & study=="Misame" & visit==vv]
  d[, feature := feature_key(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]         # smallest-p row wins per feature
  d[, .(feature, est, pval, se = se_of(cil, ciu))]
}
# Same shape as get_milk(), but for one blood dataset (dd) at one visit (vv).
get_blood <- function(dd, vv){
  d <- bloodC[measure=="ATE" & dataset==dd & visit==vv]
  d[, feature := feature_key(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]         # smallest-p row wins per feature
  d[, .(feature, est, pval, se = se_of(cil, ciu))]
}

mk_mz <- mzrt_milk(misame_only=FALSE); pl_mz <- mzrt_rlc("ProcessedDataMISAME3_plasma.csv"); vm_mz <- mzrt_vams()

# unique-mass reciprocal subset: keep candidate pairs where feature_A and feature_B
# each have exactly ONE partner within 25 ppm (+mode) — drops ambiguous collisions.
unique_match <- function(A, B){
  h <- match_mzrt(A, B, same_mode=TRUE, use_rt=FALSE)
  if(!nrow(h)) return(h)
  h[, n_partners_A := .N, by=feature_A]                 # how many B masses this A matched
  h[, n_partners_B := .N, by=feature_B]                 # how many A masses this B matched
  h[n_partners_A==1 & n_partners_B==1, .(feature_A, feature_B)]
}

# Concordance metrics for one matched skeleton (skel = feature_A<->feature_B pairs),
# joining each side's effects (fA, fB) back on.
mets <- function(skel, fA, fB){
  dt <- merge(skel, fA[,.(feature_A=feature, est_A=est, pval_A=pval, se_A=se)], by="feature_A")
  dt <- merge(dt,  fB[,.(feature_B=feature, est_B=est, pval_B=pval, se_B=se)], by="feature_B")
  dt <- dt[is.finite(est_A) & is.finite(est_B)]

  # "Anchor" pairs = those already significant on side A; concordance is more
  # meaningful among features BEP actually moved than among the noisy remainder.
  anc <- dt[pval_A < 0.05]

  # Inverse-variance weights for the effect correlation (tighter CIs count more).
  w <- 1/(dt$se_A^2 + dt$se_B^2); w[!is.finite(w)] <- NA

  data.table(n_pairs=nrow(dt),
             pct_conc=round(100*mean(sign(dt$est_A)==sign(dt$est_B))),
             n_anchor=nrow(anc),
             pct_conc_anchor=if(nrow(anc)>=5) round(100*mean(sign(anc$est_A)==sign(anc$est_B))) else NA_real_,
             w_r=round(weighted_cor(dt$est_A, dt$est_B, w), 3))
}

out <- list()
add <- function(pair, variant, skel, fA, fB) out[[length(out)+1]] <<- cbind(pair=pair, variant=variant, mets(skel, fA, fB))

# (1) milk <-> maternal plasma  (V1<->V1, same platform)
mk <- get_milk("1-2 mo."); pl <- get_blood("MaternalPlasma","pn12")
A <- mk_mz[feature %in% mk$feature]; B <- pl_mz[feature %in% pl$feature]
add("Milk <-> maternal plasma (same platform)", "mass only (25 ppm)",   best_match(A,B,use_rt=FALSE),               mk, pl)
add("Milk <-> maternal plasma (same platform)", "mass + RT (0.20 min)", best_match(A,B,rt_tol=CC_RT_V1V1,use_rt=TRUE), mk, pl)

# (2) milk <-> infant  (V1<->V3, cross platform)
inf <- get_blood("VamsPostnatalInfant","pn12")
A <- mk_mz[feature %in% mk$feature]; B <- vm_mz[feature %in% inf$feature]
add("Milk <-> infant blood (cross platform)", "mass only (25 ppm)",       best_match(A,B,use_rt=FALSE), mk, inf)
add("Milk <-> infant blood (cross platform)", "unique mass (1:1 within 25 ppm)", unique_match(A,B),     mk, inf)

# (2) plasma <-> infant  (V1<->V3, cross platform)
A <- pl_mz[feature %in% pl$feature]; B <- vm_mz[feature %in% inf$feature]
add("Maternal plasma <-> infant blood (cross platform)", "mass only (25 ppm)",       best_match(A,B,use_rt=FALSE), pl, inf)
add("Maternal plasma <-> infant blood (cross platform)", "unique mass (1:1 within 25 ppm)", unique_match(A,B),     pl, inf)

res <- rbindlist(out)
fwrite(res, paste0(root,"results/cross_compartment_matching_sensitivity.csv"))
cat("\n=== matching sensitivity (combined-arms, any-postnatal-BEP) ===\n"); print(res)
cat("\nSaved results/cross_compartment_matching_sensitivity.csv\n")
