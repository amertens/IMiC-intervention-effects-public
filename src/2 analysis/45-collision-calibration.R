# =============================================================================
# 45-collision-calibration.R
#
# Calibrates the 25 ppm cross-compartment mass-matching against chance, which
# the report (Methods item 5) flags as missing: "a one-line estimate of the mean
# number of background features within 25 ppm of a random query mass would let
# every by-mass count be calibrated."
#
# Uses ONLY the m/z-RT feature catalogues already on disk (no ATE results, no raw
# spectra). Matching is m/z within 25 ppm + same ionization mode (RT off), i.e.
# exactly the first-pass matcher used in scripts 14/19/20 (_blood_helpers.R).
#
#  Part A  Within-catalogue collision density: for each feature, how many OTHER
#          same-mode features fall within +/-25 ppm. This is the literal
#          "mean background features within 25 ppm of a random query mass".
#  Part B  Cross-catalogue rate: mean number of platform-B features within 25 ppm
#          of a random platform-A feature (same mode), for the pairs used in the
#          section-3 concordance table.
#  Part C  Permutation null for a jointly-significant by-mass count. Using the
#          significant-set SIZES already reported (NOT altered here) as inputs,
#          how many 25 ppm both-matches arise when the significant labels are
#          randomised? Demonstrated for milk<->maternal plasma (report: 936 milk
#          FDR-sig, 593 plasma FDR-sig, observed 50 matched-both).
#
# Out: results/collision_calibration.csv
# =============================================================================

suppressMessages({library(data.table); library(here)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))   # match_mzrt(), CC_PPM=25

set.seed(20240705)

cats <- list(
  Milk_V1           = tryCatch(mzrt_milk(misame_only = TRUE),                 error = function(e) NULL),
  MaternalPlasma_V1 = tryCatch(mzrt_rlc("ProcessedDataMISAME3_plasma.csv"),    error = function(e) NULL),
  PrenatalVAMS_V1   = tryCatch(mzrt_rlc("ProcessedDataMISAME3_VAMS.csv"),      error = function(e) NULL),
  PostnatalVAMS_V3  = tryCatch(mzrt_vams(),                                    error = function(e) NULL))
cats <- Filter(Negate(is.null), cats)
for (nm in names(cats)) cats[[nm]] <- cats[[nm]][!is.na(mz) & !is.na(mode)]

# -- Part A: within-catalogue neighbours within 25 ppm (same mode) ------------
# Match a catalogue against itself, drop self-pairs, then count how many OTHER
# same-mode features fall within 25 ppm of each feature. Features with no neighbour
# are absent from the pair table, so left-join back and fill their count with 0.
partA <- rbindlist(lapply(names(cats), function(nm) {
  A <- cats[[nm]]
  pairs_25ppm <- match_mzrt(A, A)             # same-mode pairs within 25 ppm (includes self)
  pairs_25ppm <- pairs_25ppm[feature_A != feature_B]
  neighbour_count <- pairs_25ppm[, .N, by = feature_A]
  per_feature <- merge(data.table(feature_A = A$feature), neighbour_count, by = "feature_A", all.x = TRUE)
  per_feature[is.na(N), N := 0L]
  data.table(platform = nm, n_features = nrow(A),
             mean_neighbours = round(mean(per_feature$N), 3),
             median_neighbours = median(per_feature$N),
             p90_neighbours = as.numeric(quantile(per_feature$N, 0.90)),
             max_neighbours = max(per_feature$N),
             frac_with_ge1 = round(mean(per_feature$N >= 1), 3))
}))

cat("\n===== Part A: within-catalogue features within +/-25 ppm (same mode) =====\n")
print(partA)

# -- Part B: cross-catalogue rate for the section-3 pairs ---------------------
pairs <- list(
  c("Milk_V1", "MaternalPlasma_V1"),   # V1<->V1 (same method)
  c("Milk_V1", "PostnatalVAMS_V3"),    # V1<->V3 (cross method) = milk -> infant blood
  c("MaternalPlasma_V1", "PostnatalVAMS_V3"))
partB <- rbindlist(lapply(pairs, function(p) {
  A <- cats[[p[1]]]; B <- cats[[p[2]]]
  if (is.null(A) || is.null(B)) return(NULL)
  h <- match_mzrt(A, B)                        # all A<->B pairs within 25 ppm, same mode
  data.table(A = p[1], B = p[2], nA = nrow(A), nB = nrow(B),
             total_pairs_25ppm = nrow(h),
             mean_B_within_25ppm_per_A = round(nrow(h) / nrow(A), 3))
}))
cat("\n===== Part B: platform-B features within 25 ppm of a random platform-A feature =====\n")
print(partB)

# -- Part C: permutation null for a jointly-significant by-mass count ----------
# Inputs are the significant-set SIZES already reported (not recomputed / not altered):
#   milk pooled FDR-sig = 936 ; maternal plasma pn12 FDR-sig = 593 ; observed both = 50.
# Null distribution of the both-significant match count: draw random "significant"
# subsets of the true sizes from each catalogue and count their 25-ppm matches. Repeating
# this n_perm times gives the number of matches expected if significance ignored mass.
perm_expected_matches <- function(A, B, nA_sig, nB_sig, n_perm = 1000) {
  nA_sig <- min(nA_sig, nrow(A)); nB_sig <- min(nB_sig, nrow(B))
  vapply(seq_len(n_perm), function(i) {
    A_rand <- A[sample.int(nrow(A), nA_sig)]
    B_rand <- B[sample.int(nrow(B), nB_sig)]
    nrow(match_mzrt(A_rand, B_rand))
  }, numeric(1))
}
partC <- NULL
if (all(c("Milk_V1", "MaternalPlasma_V1") %in% names(cats))) {
  # Use the SAME sig-set sizes that entered the matcher (per-feature-deduped FDR-sig
  # counts), read from the ATE results -- NOT the report's headline milk count, which
  # is not per-feature-deduped. Observed both-sig is read from the script-14 match file.
  count_sig <- function(rds, dataset_val = NULL, visits = NULL, study_val = NULL) {
    d <- as.data.table(readRDS(paste0(root, "results/", rds)))[measure == "ATE"]
    if (!is.null(study_val))   d <- d[study == study_val]
    if (!is.null(dataset_val)) d <- d[dataset == dataset_val]
    if (!is.null(visits))      d <- d[visit %in% visits]
    d[, feature := toupper(biomarker)]; d <- d[order(pval)][!duplicated(feature)]
    sum(d$sigFDR == 1, na.rm = TRUE)
  }
  MILK_SIG   <- count_sig("adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS", study_val = "Misame")
  PLASMA_SIG <- count_sig("blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS", dataset_val = "MaternalPlasma", visits = "pn12")
  mpcsv <- paste0(root, "results/cross_compartment_metab_MilkVsMaternalPlasma_adjusted.csv")
  OBSERVED_BOTH <- if (file.exists(mpcsv)) { dd <- fread(mpcsv); sum(dd$sigFDR_A == 1 & dd$sigFDR_B == 1, na.rm = TRUE) } else NA_integer_
  cat(sprintf("  milk-plasma inputs: milk FDR-sig=%d, plasma FDR-sig=%d, observed both=%d\n", MILK_SIG, PLASMA_SIG, OBSERVED_BOTH))
  null_counts <- perm_expected_matches(cats$Milk_V1, cats$MaternalPlasma_V1, MILK_SIG, PLASMA_SIG, 1000)
  emp_p <- (1 + sum(null_counts >= OBSERVED_BOTH)) / (length(null_counts) + 1)
  partC <- data.table(
    pair = "Milk_V1 <-> MaternalPlasma_V1",
    milk_sig_in = MILK_SIG, plasma_sig_in = PLASMA_SIG, observed_both_reported = OBSERVED_BOTH,
    expected_by_chance_mean = round(mean(null_counts), 2),
    expected_by_chance_sd = round(sd(null_counts), 2),
    null_p95 = as.numeric(quantile(null_counts, 0.95)),
    empirical_p_obs_ge = signif(emp_p, 3))
  cat("\n===== Part C: permutation null for milk<->plasma jointly-significant matched pairs =====\n")
  cat(sprintf("  Inputs (from report, not recomputed): milk FDR-sig=%d, plasma FDR-sig=%d, observed matched-both=%d\n",
              MILK_SIG, PLASMA_SIG, OBSERVED_BOTH))
  print(partC)
  cat(sprintf("\n  Interpretation: if significance were unrelated to mass, ~%.1f of the 25 ppm matches\n  would be chance collisions (95th pct %.0f); observed = %d, empirical p(obs >= chance) = %.3g.\n",
              mean(null_counts), quantile(null_counts, 0.95), OBSERVED_BOTH, emp_p))
}

# -- write --------------------------------------------------------------------
out <- rbindlist(list(
  partA[, .(part = "A_within_catalogue", item = platform, metric = "mean_features_within_25ppm", value = mean_neighbours)],
  partA[, .(part = "A_within_catalogue", item = platform, metric = "frac_features_with_ge1_neighbour", value = frac_with_ge1)],
  partB[, .(part = "B_cross_catalogue", item = paste(A, B, sep = " -> "), metric = "mean_B_within_25ppm_per_A", value = mean_B_within_25ppm_per_A)],
  partB[, .(part = "B_cross_catalogue", item = paste(A, B, sep = " -> "), metric = "total_pairs_within_25ppm", value = total_pairs_25ppm)]
), use.names = TRUE)
if (!is.null(partC))
  out <- rbindlist(list(out,
    data.table(part = "C_permutation_null", item = partC$pair, metric = "expected_both_by_chance_mean", value = partC$expected_by_chance_mean),
    data.table(part = "C_permutation_null", item = partC$pair, metric = "observed_both_reported",       value = partC$observed_both_reported),
    data.table(part = "C_permutation_null", item = partC$pair, metric = "empirical_p_obs_ge_chance",    value = partC$empirical_p_obs_ge)),
    use.names = TRUE)
fwrite(out, paste0(root, "results/collision_calibration.csv"))
cat("\nSaved: results/collision_calibration.csv\n")
