# =============================================================================
# _blood_helpers.R  —  shared helpers for the blood / cross-compartment / mummichog pipeline
#
# WHY THIS FILE EXISTS
#   The m/z–RT extractors, the cross-compartment matcher, the threshold-free
#   concordance methods, and the bioTMLE result-tidier were originally copy-pasted
#   into several numbered scripts. They drifted out of sync (different names, and a
#   foverlaps column swap that produced silent zero-matches in one copy). Everything
#   shared now lives here, once, so the numbered scripts stay short and consistent.
#
# USAGE
#   source(paste0(here::here(), "/src/2 analysis/_blood_helpers.R"))
#   (Source 0-config.R first if you need extract_blood_results(), which reuses
#    extract_res() and ci_to_pvalue() from the config.)
#
# CONTENTS
#   1. Tolerances / constants
#   2. m/z–RT extractors            (feature id -> mz, rt, ionization mode)
#   3. Cross-compartment matcher    (match_mzrt = all pairs; best_match = greedy 1:1)
#   4. Threshold-free concordance   (signed_stat, gsea_es, perm_gsea_p, rrho_directional)
#   5. bioTMLE result tidier        (extract_blood_results -> long, FDR per visit×dataset)
# =============================================================================

suppressMessages({library(data.table)})
.bh_datadir <- function() paste0(here::here(), "/data/additional datasets/")

# -- 1. tolerances / constants ------------------------------------------------
# Cross-compartment feature matching (Kim x Trenton meeting, 2026-06-24).
# RT is method-aware: V1<->V1 pairs (milk / plasma / prenatal-VAMS, same LCV1 method)
# align ~linearly so use a tight window; V1<->V3 pairs (anything <-> postnatal VAMS,
# LCB3) drift non-linearly so use a loose window and lean on m/z.
CC_PPM         <- 25     # accurate-mass tolerance, +/- ppm
CC_RT_V1V1     <- 0.20   # tight RT window (same method)
CC_RT_V1V3     <- 0.50   # loose RT window (cross method)
CC_EARLY_RT    <- 0.10   # Kim false-positive trim: drop a pair if it elutes before this...
CC_EARLY_SHIFT <- 0.10   # ...AND its two RTs differ by more than this

# Mummichog (matches the milk pipeline / Trenton's recipe). Note: MUM_PPM is the
# database annotation tolerance and is a DIFFERENT quantity from CC_PPM above.
MUM_PPM    <- 10
MUM_NET    <- "human_mfn"
MUM_CUTOFF <- 0.05

# -- 2. m/z–RT extractors -----------------------------------------------------
# All return data.table(feature [UPPERCASE id], mz, rt, mode).
ion_of <- function(x) ifelse(grepl("_POS", toupper(x)), "positive",
                      ifelse(grepl("_NEG", toupper(x)), "negative", NA_character_))

# rLC catalogues — maternal plasma / prenatal VAMS — from the processed Sapient tables.
mzrt_rlc <- function(file) {
  d <- fread(paste0(.bh_datadir(), file), select = c("MZ", "RT", "Metabolite_Feature_Label"))
  data.table(feature = toupper(d$Metabolite_Feature_Label), mz = d$MZ, rt = d$RT,
             mode = ion_of(d$Metabolite_Feature_Label))[!is.na(mz)]
}
# postnatal VAMS — maternal + infant share this one V3 `vam_` catalogue.
# drop_v3_lipids=TRUE removes the di/triglyceride & cholesterol-ester features that the
# V3 method captures but V1 does not (Kim's V1<->V3 false-positive trim). Only 37 of
# 38,761 features carry such a class label, and none are in any current matched pair,
# so this has no effect on present results — provided for completeness.
mzrt_vams <- function(drop_v3_lipids = FALSE) {
  d <- fread(paste0(.bh_datadir(), "metabolite_description_vam_with_global_id.csv"))
  if (drop_v3_lipids)
    d <- d[!grepl("triglycer|diglycer|diacylglyc|triacylglyc|cholesteryl|cholesterol ester|acylglycerol",
                  d[["Compound Class"]], ignore.case = TRUE)]
  data.table(feature = toupper(d$feature_label), mz = d$mz, rt = d$rt_minute,
             mode = tolower(d$ionization_mode))[!is.na(mz) & !is.na(rt)]
}
# milk — IMiC cross-study alignment key. misame_only=TRUE keeps just the MISAME-3 ids
# (for matching against MISAME-3 blood); FALSE also includes CHILD/ELICIT/VITAL.
mzrt_milk <- function(misame_only = TRUE) {
  a <- fread(paste0(.bh_datadir(), "IMiC_alignment.csv"))
  d <- data.table(feature = toupper(a$mtb_id_MISAME3), mz = a$mz_MISAME3, rt = a$rt_MISAME3)
  if (!misame_only)
    d <- rbind(d, data.table(feature = toupper(a$mtb_id_CHILD_ELICIT_VITAL),
                             mz = a$mz_CHILD_ELICIT_VITAL, rt = a$rt_CHILD_ELICIT_VITAL))
  d <- d[!is.na(mz) & !is.na(rt) & feature != ""]
  d[, mode := ion_of(feature)]; unique(d)
}

# -- 3. cross-compartment matcher ---------------------------------------------
# PPM match = Kim's formula: |m_A − m_B| / m * 1e6 ≤ CC_PPM (±25 ppm), implemented as
# the window m_A*(1 ± CC_PPM/1e6). RT is an OPTIONAL confidence filter — per the
# 2026-06-24 meeting (Trenton: "wouldn't worry about retention time right now"),
# the FIRST-PASS match is m/z + ionization mode only (use_rt = FALSE). Pass
# use_rt = TRUE to add the method-aware RT window + Kim's early-elution trim for a
# higher-confidence list.
#
# foverlaps note: foverlaps(x = B, y = A) returns A's columns UNPREFIXED and B's with
# an `i.` prefix. So `feature` is A's, `i.feature` is B's. (Swapping these silently
# produced zero matches once — keep it here, once.)
.mz_join <- function(A, B, rt_tol, same_mode = TRUE, use_rt = FALSE) {
  A <- copy(A); B <- copy(B)
  A[, `:=`(mz_lo = mz*(1 - CC_PPM/1e6), mz_hi = mz*(1 + CC_PPM/1e6))]
  B[, `:=`(mz_lo = mz, mz_hi = mz)]; setkey(A, mz_lo, mz_hi); setkey(B, mz_lo, mz_hi)
  h <- foverlaps(B, A, type = "within", nomatch = 0L)
  if (same_mode && all(c("mode", "i.mode") %in% names(h)))
    h <- h[!is.na(mode) & !is.na(i.mode) & mode == i.mode]
  if (use_rt) {                                                          # optional RT confidence layer
    h <- h[abs(rt - i.rt) <= rt_tol]
    h <- h[!(pmin(rt, i.rt) < CC_EARLY_RT & abs(rt - i.rt) > CC_EARLY_SHIFT)]   # Kim FP trim
  }
  h
}
# every A<->B feature pair within tolerance
match_mzrt <- function(A, B, rt_tol = CC_RT_V1V3, same_mode = TRUE, use_rt = FALSE) {
  h <- .mz_join(A, B, rt_tol, same_mode, use_rt)
  if (!nrow(h)) return(data.table(feature_A = character(), feature_B = character()))
  unique(h[, .(feature_A = feature, feature_B = i.feature)])
}
# greedy 1:1 — each A feature -> its nearest-ppm B, then each B used once
best_match <- function(A, B, rt_tol = CC_RT_V1V3, same_mode = TRUE, use_rt = FALSE) {
  h <- .mz_join(A, B, rt_tol, same_mode, use_rt)
  if (!nrow(h)) return(data.table(feature_A = character(), feature_B = character()))
  h[, ppm := abs(mz - i.mz)/i.mz*1e6]; setorder(h, ppm)
  h <- h[!duplicated(feature)][!duplicated(i.feature)]
  h[, .(feature_A = feature, feature_B = i.feature)]
}

# -- 4. threshold-free concordance methods ------------------------------------
# signed significance: positive = up-regulated, magnitude = -log10 p
signed_stat <- function(pval, est) (-log10(pmax(pval, 1e-300))) * sign(est)

# GSEA running-sum enrichment score of `in_set` within a ranked `stat` vector.
# Walk features from most to least positive `stat`; a running sum steps UP at each
# set member (weighted by |stat|) and DOWN at each non-member. The score is the
# largest deviation of that walk from zero.
gsea_es <- function(stat, in_set) {
  rank_order <- order(stat, decreasing = TRUE)
  is_member  <- in_set[rank_order]                 # is each ranked feature in the set?
  weight     <- abs(stat[rank_order])              # step size = magnitude of the ranking stat
  n_members  <- sum(is_member)
  if (n_members == 0 || n_members == length(is_member)) return(0)  # no enrichment if all/none in set
  hit_fraction  <- cumsum(ifelse(is_member, weight, 0)) / sum(weight[is_member])
  miss_fraction <- cumsum(ifelse(!is_member, 1, 0)) / (length(is_member) - n_members)
  running_sum   <- hit_fraction - miss_fraction
  running_sum[which.max(abs(running_sum))]         # peak deviation = enrichment score
}
# Permutation p-value (+ normalized ES) for the enrichment of `in_set`: compare the
# observed score against scores of random sets of the same size. Returns list(es, nes, p).
perm_gsea_p <- function(stat, in_set, P = 1000) {
  observed <- gsea_es(stat, in_set)
  set_size <- sum(in_set)
  n_total  <- length(in_set)
  if (set_size < 3) return(list(es = observed, nes = NA_real_, p = NA_real_))  # too small to permute
  null <- vapply(seq_len(P), function(i) {
    random_set <- logical(n_total)
    random_set[sample.int(n_total, set_size)] <- TRUE
    gsea_es(stat, random_set)
  }, numeric(1))
  same_sign <- null[sign(null) == sign(observed)]  # null scores pointing the observed direction
  list(es  = observed,
       nes = if (length(same_sign)) observed / mean(abs(same_sign)) else NA_real_,
       p   = (1 + sum(abs(null) >= abs(observed))) / (P + 1))   # +1 smoothing avoids p = 0
}
# RRHO directional overlap of two signed-stat vectors over a COMMON item set.
# Returns c(concordant, discordant) = peak -log10 hypergeometric p in the
# up-up/down-down vs up-down/down-up corners (top half only, to avoid the
# trivial whole-set inflation).
rrho_directional <- function(sA, sB, steps = 40) {
  n <- length(sA)
  if (n < 20) return(c(concordant = NA_real_, discordant = NA_real_))
  cuts <- round(seq(1, n/2, length.out = steps))    # thresholds scanned in each list (top half only)
  ranked <- function(s, up) rank(if (up) -s else s) # rank up-most (up=TRUE) or down-most first
  # Peak overlap significance in one RRHO corner: over every pair of top-a / top-b
  # thresholds, score the hypergeometric enrichment of shared features and keep the max
  # (as -log10 p).
  corner <- function(rank_A, rank_B) {
    peak <- 0
    for (a in cuts) {
      in_top_A <- rank_A <= a
      for (b in cuts) {
        overlap <- sum(in_top_A & rank_B <= b)
        neg_log10_p <- -phyper(overlap - 1, b, n - b, a, lower.tail = FALSE, log.p = TRUE) / log(10)
        peak <- max(peak, neg_log10_p)
      }
    }
    peak
  }
  up_up     <- corner(ranked(sA, TRUE),  ranked(sB, TRUE))    # both up-regulated
  down_down <- corner(ranked(sA, FALSE), ranked(sB, FALSE))   # both down-regulated
  up_down   <- corner(ranked(sA, TRUE),  ranked(sB, FALSE))   # A up while B down
  down_up   <- corner(ranked(sA, FALSE), ranked(sB, TRUE))    # A down while B up
  c(concordant = max(up_up, down_down), discordant = max(up_down, down_up))
}
# inverse-variance-weighted Pearson r (w = 1/se^2); falls back to unweighted
weighted_cor <- function(x, y, w) {
  if (sum(is.finite(w)) <= 10) return(cor(x, y))   # too few usable weights -> plain Pearson r
  weighted_mean <- function(v) sum(w * v, na.rm = TRUE) / sum(w, na.rm = TRUE)
  mx <- weighted_mean(x); my <- weighted_mean(y)
  covariance <- sum(w * (x - mx) * (y - my), na.rm = TRUE)
  var_x      <- sum(w * (x - mx)^2, na.rm = TRUE)
  var_y      <- sum(w * (y - my)^2, na.rm = TRUE)
  covariance / sqrt(var_x * var_y)
}
# Deming (orthogonal) regression slope, assuming equal error variances
deming_slope <- function(x, y) {
  sxx <- var(x); syy <- var(y); sxy <- cov(x, y)
  (syy - sxx + sqrt((syy - sxx)^2 + 4*sxy^2)) / (2*sxy)
}

# -- 5. bioTMLE result tidier -------------------------------------------------
# Raw per-compartment bioTMLE output (named list) -> long data.frame with
# Benjamini-Hochberg FDR applied PER visit x dataset (group_by studytime, measure).
# Reuses extract_res() and ci_to_pvalue() from 0-config.R.
extract_blood_results <- function(res_list, compartment_label) {
  if (is.null(res_list)) return(NULL)
  ok <- !vapply(res_list$res, inherits, logical(1), what = "try-error")
  if (!any(ok)) return(NULL)
  res <- data.table::rbindlist(Map(extract_res, res_list$res[ok]), idcol = "studytime") |> as.data.frame()
  res$pval <- ci_to_pvalue(cil = res$cil, ciu = res$ciu)
  res <- res |>
    dplyr::group_by(studytime, measure) |> dplyr::mutate(pval_adj = p.adjust(pval, method = "BH")) |>
    dplyr::group_by(measure)            |> dplyr::mutate(pval_adj_global = p.adjust(pval, method = "BH")) |>
    dplyr::ungroup()
  res$sig <- 1 * (res$pval < 0.05); res$sigFDR <- 1 * (res$pval_adj < 0.05)
  res$dataset <- compartment_label; res$visit <- sub("^.*-", "", res$studytime)
  res$compartment <- compartment_label
  res
}
