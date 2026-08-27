# =============================================================================
# 45b-collision-direction-null.R   (companion to 45-collision-calibration.R)
#
# Two additions the count-only calibration in script 45 did not cover:
#
#  Part D  DIRECTION-AGREEMENT null. Among the jointly-FDR-significant matched
#          pairs, is the same-direction fraction (the report's 78% / 100%) above
#          what marginal skew alone predicts? Both milk (~38% up) and plasma
#          (~45% up) lean toward decreases, and infant VAMS is ~all up, so a high
#          same-direction fraction can arise with no paired biology. Null: hold
#          the matched pairs fixed, permute one side's signs (preserving its
#          up/down counts), recompute the concordant fraction. Also reported as
#          the analytic independence value p_A*p_B + (1-p_A)*(1-p_B).
#
#  Part E  COUNT null for the two cross-method pairs (milk->infant, plasma->
#          infant), the same permutation as script 45 Part C but for the V1<->V3
#          pairs, using the FDR-sig set sizes read from the match files.
#
# Reads the existing script-14 match outputs (matched pairs + est + sigFDR +
# both_sig) and the m/z-RT catalogues. Alters no published value; reads signs and
# counts only.
#
# Out: results/collision_calibration_direction.csv
# =============================================================================

suppressMessages({library(data.table); library(here)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))
set.seed(20240705)

resdir <- paste0(root, "results/")
pairs <- list(
  list(name = "Milk -> MaternalPlasma (V1-V1, by mass)",  csv = "cross_compartment_metab_MilkVsMaternalPlasma_adjusted.csv",  mass = TRUE,  catA = "milk", catB = "plasma"),
  list(name = "Milk -> InfantVAMS (V1-V3, by mass)",      csv = "cross_compartment_metab_MilkVsInfantVAMS_adjusted.csv",      mass = TRUE,  catA = "milk", catB = "infant"),
  list(name = "MaternalPlasma -> InfantVAMS (V1-V3, mass)", csv = "cross_compartment_metab_MaternalPlasmaVsInfantVAMS_adjusted.csv", mass = TRUE, catA = "plasma", catB = "infant"),
  list(name = "MaternalVAMS -> InfantVAMS (exact identifier)", csv = "cross_compartment_metab_MaternalVAMSvsInfantVAMS_adjusted.csv", mass = FALSE, catA = NA, catB = NA))

cat_of <- list(
  milk   = function() mzrt_milk(misame_only = TRUE),
  plasma = function() mzrt_rlc("ProcessedDataMISAME3_plasma.csv"),
  infant = function() mzrt_vams())
cat_cache <- new.env()
getcat <- function(k) { if (is.null(cat_cache[[k]])) cat_cache[[k]] <- cat_of[[k]]()[!is.na(mz) & !is.na(mode)]; cat_cache[[k]] }

as_bool <- function(x) x %in% c(TRUE, "TRUE", "True", "true", 1, "1")

# -- Part D: direction-agreement null -----------------------------------------
cat("\n===== Part D: direction-agreement among jointly-significant matched pairs =====\n")
D <- rbindlist(lapply(pairs, function(p) {
  f <- paste0(resdir, p$csv); if (!file.exists(f)) return(NULL)
  d  <- fread(f)
  both_sig_pairs <- d[as_bool(both_sig) & is.finite(est_A) & is.finite(est_B)]
  n  <- nrow(both_sig_pairs); if (n == 0) return(data.table(pair = p$name, n_both_sig = 0L))
  signA <- sign(both_sig_pairs$est_A); signB <- sign(both_sig_pairs$est_B)
  obs_same_dir <- mean(signA == signB)                    # observed concordant fraction
  # marginal up-rate on each side; if directions were independent, the concordant
  # fraction expected from skew alone is p_A*p_B + (1-p_A)*(1-p_B).
  uprate_A <- mean(both_sig_pairs$est_A > 0); uprate_B <- mean(both_sig_pairs$est_B > 0)
  analytic <- uprate_A*uprate_B + (1-uprate_A)*(1-uprate_B)
  # permutation null: keep side A fixed, reshuffle side B's signs (preserves its up/down
  # counts) and recompute concordance 5000 times.
  null <- vapply(1:5000, function(i) mean(signA == sample(signB)), numeric(1))
  data.table(pair = p$name, n_both_sig = n,
             obs_same_dir = round(obs_same_dir, 3),
             uprate_A = round(uprate_A, 3), uprate_B = round(uprate_B, 3),
             null_same_dir_analytic = round(analytic, 3),
             null_same_dir_perm = round(mean(null), 3),
             emp_p_obs_ge = signif((1 + sum(null >= obs_same_dir - 1e-9)) / (length(null) + 1), 3))
}))
print(D)

# -- Part E: count null for the two cross-method (V1<->V3) pairs ---------------
cat("\n===== Part E: count null (chance 25-ppm both-significant matches), cross-method pairs =====\n")
# Same count-null idea as script 45 Part C: draw random significant subsets of the true
# sizes and count their 25-ppm matches, so `obs` can be compared against chance collisions.
count_null <- function(A, B, nA_sig, nB_sig, obs, n_perm = 1000) {
  nA_sig <- min(nA_sig, nrow(A)); nB_sig <- min(nB_sig, nrow(B))
  null <- vapply(seq_len(n_perm), function(i)
    nrow(match_mzrt(A[sample.int(nrow(A), nA_sig)], B[sample.int(nrow(B), nB_sig)])), numeric(1))
  list(mean = mean(null), sd = sd(null), p95 = quantile(null, 0.95),
       emp_p = (1 + sum(null >= obs)) / (n_perm + 1))
}
# TRUE marginal FDR-sig set sizes per compartment. Deriving these from the matched-pairs
# CSV would exclude every FDR-sig feature with no 25-ppm partner and bias the null low,
# so read the totals straight from the ATE results.
count_sig <- function(rds, dataset_val = NULL, visits = NULL, study_val = NULL) {
  d <- as.data.table(readRDS(paste0(resdir, rds)))[measure == "ATE"]
  if (!is.null(study_val))   d <- d[study == study_val]
  if (!is.null(dataset_val)) d <- d[dataset == dataset_val]
  if (!is.null(visits))      d <- d[visit %in% visits]
  d[, feature := toupper(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]                # one row per feature
  sum(d$sigFDR == 1, na.rm = TRUE)
}
sig_tot <- list(
  milk   = count_sig("adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS", study_val = "Misame"),
  plasma = count_sig("blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS", dataset_val = "MaternalPlasma", visits = "pn12"),
  infant = count_sig("blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS", dataset_val = "VamsPostnatalInfant", visits = c("pn12","pn34","pn56")))
cat(sprintf("  marginal FDR-sig totals: milk=%d, plasma=%d, infant=%d\n", sig_tot$milk, sig_tot$plasma, sig_tot$infant))
E <- rbindlist(lapply(pairs[2:3], function(p) {
  f <- paste0(resdir, p$csv); if (!file.exists(f)) return(NULL)
  d <- fread(f)
  nA_sig <- sig_tot[[p$catA]]; nB_sig <- sig_tot[[p$catB]]
  obs    <- nrow(d[as_bool(both_sig)])
  cn <- count_null(getcat(p$catA), getcat(p$catB), nA_sig, nB_sig, obs)
  data.table(pair = p$name, sigA_total = nA_sig, sigB_total = nB_sig,
             observed_both = obs, expected_by_chance = round(cn$mean, 2),
             chance_sd = round(cn$sd, 2), chance_p95 = as.numeric(round(cn$p95, 1)),
             emp_p_obs_ge = signif(cn$emp_p, 3))
}))
print(E)

fwrite(rbindlist(list(
  D[, .(part = "D_direction_null", pair, n = n_both_sig, metric = "observed_same_direction", value = obs_same_dir)],
  D[, .(part = "D_direction_null", pair, n = n_both_sig, metric = "null_same_direction_perm", value = null_same_dir_perm)],
  D[, .(part = "D_direction_null", pair, n = n_both_sig, metric = "emp_p_obs_ge",             value = emp_p_obs_ge)],
  E[, .(part = "E_count_null_crossmethod", pair, n = observed_both, metric = "expected_both_by_chance", value = expected_by_chance)],
  E[, .(part = "E_count_null_crossmethod", pair, n = observed_both, metric = "emp_p_obs_ge",           value = emp_p_obs_ge)]
), use.names = TRUE), paste0(resdir, "collision_calibration_direction.csv"))
cat("\nSaved: results/collision_calibration_direction.csv\n")
