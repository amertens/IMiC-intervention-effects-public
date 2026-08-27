# =============================================================================
# 40-dyadic-blood-distance.R
#
# "Beta-diversity"-style DYADIC compositional similarity of MOTHER vs INFANT blood,
# tested by BEP arm (April's idea, 2026-06-30 IMiC meeting; Andrew action item).
#
# RATIONALE: rather than tracking individual transferred features, ask whether BEP
# makes a mother's and her OWN infant's blood metabolome more ALIKE. Postnatal VAMS is
# the only clean place for this: mother & infant share the SAME V3 `vam_` catalogue
# (identity-anchored -- no cross-platform m/z matching), and the parsed `idBiospe` is
# the FAMILY key shared by a mother and her infant (the prep joins BOTH the mother and
# her infant rows to the treatment arm by idBiospe, which only works because it is
# family-level) -> genuine mother-infant dyads.
#
# Per dyad (a mother row + an infant row at the SAME visit and SAME idBiospe) we
# collapse to ONE similarity over the shared feature set, then test whether it differs
# by arm. HYPOTHESIS: dyad similarity is HIGHER (distance LOWER) under postnatal BEP,
# strongest for the postnatal-only contrast.
#
# Metrics per dyad:
#   euclid   Euclidean distance on the (already z-scored) intensities  (lower = closer)
#   cor_sim  Spearman rho between the two profiles                       (higher = closer)
# computed over (a) ALL shared features and (b) the BEP-affected subset (infant FDR-sig
# ATE features) -- Andrew's "do it within the components BEP actually moved".
#
# NOTE on timepoints: maternal postnatal VAMS exists only at tri3/pn56, so dyads can
# only form where both sides were sampled (mostly pn56) -- which is exactly the
# identity-anchored transfer timepoint. The by-visit table shows where dyads exist.
#
# Out: results/dyadic_blood_distance.csv          (one row per dyad x visit x feature-set)
#      results/dyadic_blood_distance_armtest.csv   (arm-effect tests per visit x feature-set)
#      figures/cross_compartment/dyadic_distance_by_arm.png
# =============================================================================

suppressMessages({library(data.table); library(ggplot2)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))
set.seed(20240630)
outfig <- paste0(root, "figures/cross_compartment"); dir.create(outfig, showWarnings = FALSE, recursive = TRUE)

if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- TRUE
.btag <- if (BLOOD_ADJUST) "adjusted_" else ""

# -- postnatal VAMS (mother + infant on one catalogue) -------------------------
mb <- readRDS(paste0(root, "data/blood/merged_blood_datasets.RDS"))
pv <- as.data.table(mb$vams_postnatal)
featcols <- grep("^vam_", names(pv), value = TRUE, ignore.case = TRUE)
if (length(featcols) < 50)
  warning("Only ", length(featcols), " vam_ feature columns in merged_blood_datasets.RDS -- ",
          "looks like a SETUP (100-feature) prep. Re-run the FULL prep for the real metric.")
cat(sprintf("postnatal VAMS rows=%d  features=%d  dyad={%s}  visits={%s}\n",
            nrow(pv), length(featcols), paste(unique(pv$dyad), collapse = ","),
            paste(sort(unique(pv$visit)), collapse = ",")))

# -- pair mother & infant by family (idBiospe) x visit -------------------------
# A dyad is one family (idBiospe) sampled at one visit; the key "idBiospe|visit"
# identifies a mother row and her infant row that belong together. Keep one row per
# (family, visit) on each side, then intersect keys to get families sampled on BOTH sides.
key  <- c("idBiospe", "visit")
moth <- pv[dyad == "mother"]; inf <- pv[dyad == "infant"]
moth <- moth[!duplicated(moth[, ..key])]
inf  <- inf [!duplicated(inf [, ..key])]
Mmat <- as.matrix(moth[, ..featcols]); rownames(Mmat) <- paste(moth$idBiospe, moth$visit, sep = "|")
Imat <- as.matrix(inf [, ..featcols]); rownames(Imat) <- paste(inf$idBiospe,  inf$visit,  sep = "|")
common <- intersect(rownames(Mmat), rownames(Imat))     # dyad keys present on both mother and infant sides
arm_lkp <- unique(moth[, .(dyad_key = paste(idBiospe, visit, sep = "|"), arm, class)])
cat(sprintf("Mother-infant dyads (idBiospe x visit on both sides): %d  (families %d)\n",
            length(common), uniqueN(sub("\\|.*", "", common))))
if (!length(common)) stop("No mother-infant dyads matched -- check idBiospe parsing / visit overlap.")

# -- BEP-affected feature subset (infant FDR-sig ATE features) -----------------
bloodS <- as.data.table(readRDS(paste0(root, "results/blood_compartment_", .btag,
            "intervention_effects_results_clean.RDS")))
bep_feats <- toupper(unique(bloodS[dataset == "VamsPostnatalInfant" & measure == "ATE" & sigFDR == 1]$biomarker))
sub_idx   <- which(toupper(featcols) %in% bep_feats)
cat(sprintf("BEP-affected (infant FDR-sig) features available in matrix: %d\n", length(sub_idx)))

# -- per-dyad metrics over a feature subset ------------------------------------
# For every dyad, collapse the mother profile x and infant profile y (over the chosen
# feature columns) to two similarity numbers: Euclidean distance (lower = closer) and
# Spearman correlation (higher = closer). `idx = NULL` means use all shared features.
dyad_metrics <- function(idx, label) {
  cols <- if (is.null(idx)) seq_along(featcols) else idx
  if (length(cols) < 3) { message("skip '", label, "': <3 features"); return(NULL) }
  out <- rbindlist(lapply(common, function(k) {
    mother_profile <- Mmat[k, cols]
    infant_profile <- Imat[k, cols]
    data.table(dyad_key = k,
               euclid  = sqrt(sum((mother_profile - infant_profile)^2)),
               cor_sim = suppressWarnings(cor(mother_profile, infant_profile, method = "spearman")))
  }))
  out[, feature_set := label]; out
}
dm <- rbind(dyad_metrics(NULL, "all_features"),
            dyad_metrics(sub_idx, "bep_affected"))
dm <- merge(dm, arm_lkp, by = "dyad_key")
dm[, c("idBiospe", "visit") := tstrsplit(dyad_key, "|", fixed = TRUE)]
fwrite(dm, paste0(root, "results/dyadic_blood_distance.csv"))

# -- arm-effect test: postnatal-BEP (class) contrast per visit x feature-set ---
# class = 0 control / 1 postnatal-BEP. Wilcoxon (rank-sum) tests whether dyad
# similarity differs by arm; reported for both the correlation and distance metrics.
armtest <- dm[, {
  wilcox_p_cor    <- tryCatch(wilcox.test(cor_sim ~ class)$p.value, error = function(e) NA_real_)
  wilcox_p_euclid <- tryCatch(wilcox.test(euclid  ~ class)$p.value, error = function(e) NA_real_)
  mean_cor_ctrl <- mean(cor_sim[class == 0], na.rm = TRUE)
  mean_cor_bep  <- mean(cor_sim[class == 1], na.rm = TRUE)
  .(n_dyads = .N, n_ctrl = sum(class == 0), n_bep = sum(class == 1),
    cor_ctrl = round(mean_cor_ctrl, 3), cor_bep = round(mean_cor_bep, 3),
    cor_diff = round(mean_cor_bep - mean_cor_ctrl, 3),
    wilcox_p_cor = signif(wilcox_p_cor, 3), wilcox_p_euclid = signif(wilcox_p_euclid, 3))
}, by = .(feature_set, visit)][order(feature_set, visit)]
fwrite(armtest, paste0(root, "results/dyadic_blood_distance_armtest.csv"))

cat("\n=== Dyadic mother-infant blood similarity: postnatal-BEP (class) contrast ===\n")
print(armtest)
cat("\n(cor_diff > 0 and small wilcox_p_cor => BEP raises mother-infant blood similarity)\n")

# -- 4-level arm model (richer view, pooled over visits) -----------------------
# Regress dyad similarity on the full 4-level arm (and visit) instead of the binary
# class split. Only include a predictor if it actually varies (>=2 levels), otherwise
# lm() would fail on a constant term.
for (fs in unique(dm$feature_set)) {
  d <- dm[feature_set == fs & is.finite(cor_sim)]
  preds <- c("arm", "visit")[c(uniqueN(d$arm) >= 2, uniqueN(d$visit) >= 2)]
  if (length(preds) && nrow(d) >= 10) {
    f <- reformulate(preds, "cor_sim")
    cat(sprintf("\n-- lm(%s) | feature_set = %s --\n", deparse(f), fs))
    print(summary(lm(f, data = d))$coefficients)
  }
}

# -- figure --------------------------------------------------------------------
pdat <- dm[is.finite(cor_sim)]
if (nrow(pdat)) {
  p <- ggplot(pdat, aes(arm, cor_sim, fill = arm)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_grid(feature_set ~ visit) +
    labs(x = NULL, y = "Mother-infant blood profile similarity (Spearman rho)",
         title = "Dyadic mother-infant blood similarity by BEP arm",
         subtitle = "postnatal VAMS, shared V3 catalogue (identity-anchored)") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggsave(p, file = paste0(outfig, "/dyadic_distance_by_arm.png"), width = 10, height = 6)
}

cat("\nSaved: results/dyadic_blood_distance.csv, results/dyadic_blood_distance_armtest.csv,\n",
    "       figures/cross_compartment/dyadic_distance_by_arm.png\n")
