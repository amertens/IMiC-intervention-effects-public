# =============================================================================
# 54-supplement-detection-status.R
#
# Ports the BEP/APSE supplement-DETECTION components of Trenton's "Compartment
# Tracking" script that were not yet in the repo pipeline. The repo already
# recomputes blood FDR per compartment x timepoint (script 12 / _blood_helpers.R),
# runs directional Mummichog (18/25), annotates putative names (23/47), and does
# 25-ppm cross-compartment mass matching (14/39). What was MISSING and is added
# here, faithful to Trenton's Compartment_Tracking.Rmd:
#
#   (4) raw APSE supplement DETECTION profile from ProcessedDataMISAME3_VAMS.csv
#       - a feature is "reliably detected" when present in >= half of the 18 APSE
#         replicates; "supplement-abundant" when its mean abundance is >= 1 SD
#         above the mean abundance ACROSS supplement features.
#   (5) hardcoded significant-count validation (stopifnot) against the counts
#       Trenton reported for the blood compartments.
#   (8) graded 5-level supplement_status for every FDR-significant feature
#       (No matched supplement feature / Not detected / Detected /
#        Reliably detected / Supplement-abundant), matched exact-id first then
#        same-ion-mode within 25 ppm.
#   (9) the putative-NAME cross-compartment link tier (features linked across
#       compartments by an identical normalized putative name), direction-only.
#  (10) the two supplement figures: the supplement-status-shaded cross-compartment
#       tracking line plot, and the stacked supplement-status bar.
#
# Cross-compartment comparison is DIRECTION-ONLY (Kim intensity constraint), as
# in the rest of the blood pipeline. MISAME-III, BEP-vs-control combined arms.
#
# Out: results/bep_supplement_detection_profile.csv        (per-feature APSE tiers)
#      results/supplement_status_fdr_features.csv           (5-level status per FDR feature)
#      results/compartment_tracking_supplement.xlsx         (all sheets)
#      figures/supplement_status_by_compartment.png         (stacked bar)
#      figures/compartment_tracking_supplement_status.png   (tracking line plot)
# =============================================================================

suppressMessages({ library(data.table); library(tidyverse); library(ggrepel) })
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))

FDR_THRESHOLD <- 0.05
MASS_TOL_PPM  <- CC_PPM          # 25 ppm, the shared cross-compartment constant

# ---------------------------------------------------------------------------
# (4) Raw APSE supplement detection profile
# ---------------------------------------------------------------------------
# ProcessedDataMISAME3_VAMS.csv carries 27 "supplement;..." columns; Trenton used
# only the 18 that match "^supplement;APSE" (two APSE_Mere donors x 9 injections).
vams_file <- paste0(root, "data/additional datasets/ProcessedDataMISAME3_VAMS.csv")
supp_cols <- grep("^supplement;APSE", names(fread(vams_file, nrows = 0)), value = TRUE)
stopifnot(length(supp_cols) == 18)          # Trenton's guard: 18 APSE replicates

supp_raw <- fread(vams_file,
                  select = c("MZ", "RT", "Metabolite_Feature_Label", supp_cols))
supp_long <- melt(supp_raw,
                  id.vars = c("MZ", "RT", "Metabolite_Feature_Label"),
                  variable.name = "supplement_sample",
                  value.name = "supplement_abundance")

supplement_profile <- supp_long[, .(
  supplement_replicates       = .N,
  detected_replicates         = sum(!is.na(supplement_abundance) & supplement_abundance > 0),
  mean_supplement_abundance   = mean(supplement_abundance, na.rm = TRUE),
  median_supplement_abundance = median(supplement_abundance, na.rm = TRUE),
  sd_supplement_abundance     = sd(supplement_abundance, na.rm = TRUE)
), by = .(MZ, RT, Metabolite_Feature_Label)]

# mean/median/sd are NaN for all-missing features -> NA
for (col in c("mean_supplement_abundance", "median_supplement_abundance", "sd_supplement_abundance"))
  supplement_profile[is.nan(get(col)), (col) := NA_real_]

supplement_profile[, `:=`(
  supplement_feature           = Metabolite_Feature_Label,
  feature                      = toupper(Metabolite_Feature_Label),
  supplement_mz                = MZ,
  supplement_rt                = RT,
  ion_mode                     = ion_of(Metabolite_Feature_Label),
  supplement_detected          = detected_replicates > 0,
  supplement_reliably_detected = detected_replicates >= ceiling(supplement_replicates / 2)
)]

# "abundant" threshold = mean + 1 SD of the per-feature mean abundance, taken
# ACROSS supplement features (Trenton's definition).
abund_threshold <- supplement_profile[
  , mean(mean_supplement_abundance, na.rm = TRUE) + sd(mean_supplement_abundance, na.rm = TRUE)]
supplement_profile[, supplement_abundant :=
  supplement_reliably_detected & mean_supplement_abundance >= abund_threshold]

fwrite(supplement_profile[, .(feature, supplement_feature, supplement_mz, supplement_rt,
                              ion_mode, supplement_replicates, detected_replicates,
                              mean_supplement_abundance, median_supplement_abundance,
                              sd_supplement_abundance, supplement_detected,
                              supplement_reliably_detected, supplement_abundant)],
       paste0(root, "results/bep_supplement_detection_profile.csv"))

cat(sprintf("supplement detection profile: %d features | reliably detected %d | abundant %d\n",
            nrow(supplement_profile),
            sum(supplement_profile$supplement_reliably_detected),
            sum(supplement_profile$supplement_abundant, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# (5) Blood significant counts + hardcoded validation
# ---------------------------------------------------------------------------
# Recompute BH FDR within compartment x timepoint on the repo-generated blood
# hand-off table (script 50), exactly as Trenton did, and gate on his reported
# counts. compartment labels differ between the two scripts; map to Trenton's.
blood_in <- fread(paste0(root, "results/blood_mummichog_input_for_trenton.csv"))
blood_in[, compartment := fcase(
  compartment == "Maternal plasma",                 "Maternal plasma",
  compartment == "Maternal blood (postnatal VAMS)", "Maternal VAMS",
  compartment == "Infant blood (postnatal VAMS)",   "Infant VAMS",
  default = compartment)]
blood_in[, fdr := p.adjust(raw_p, "BH"), by = .(compartment, timepoint)]

reported_counts <- data.table(
  compartment = c("Maternal plasma", "Maternal plasma", "Maternal plasma",
                  "Maternal VAMS", "Maternal VAMS",
                  "Infant VAMS", "Infant VAMS", "Infant VAMS", "Infant VAMS"),
  timepoint   = c("incl", "tri3", "pn12", "tri3", "pn56",
                  "acco", "pn12", "pn34", "pn56"),
  reported_n  = c(17L, 10L, 593L, 0L, 18L, 0L, 3L, 8L, 18L))

calc_counts <- blood_in[, .(calc_n = sum(fdr < FDR_THRESHOLD, na.rm = TRUE)),
                        by = .(compartment, timepoint)]
count_validation <- merge(reported_counts, calc_counts,
                          by = c("compartment", "timepoint"), all.x = TRUE)
count_validation[is.na(calc_n), calc_n := 0L]
count_validation[, matches_report := calc_n == reported_n]
print(count_validation[order(compartment, timepoint)])
stopifnot(all(count_validation$matches_report))      # Trenton's hard validation gate
cat("blood significant-count validation: all 9 compartment x timepoint counts match Trenton\n")

# ---------------------------------------------------------------------------
# Assemble FDR-significant features across compartments (blood + milk)
# ---------------------------------------------------------------------------
# Blood: the compartment x timepoint FDR-significant features, with mz / ion mode
# / direction from the hand-off table and a putative name (where annotated).
blood_ann <- fread(paste0(root, "results/fdr_sig_putative_annotation.csv"))
blood_ann[, dataset := fcase(
  dataset == "MaternalPlasma",        "Maternal plasma",
  dataset == "VamsPostnatalMaternal", "Maternal VAMS",
  dataset == "VamsPostnatalInfant",   "Infant VAMS",
  default = dataset)]
blood_ann[, feature := toupper(feature)]
blood_ann_key <- unique(blood_ann[, .(compartment = dataset, timepoint = visit,
                                       feature = feature, putative_name = best_annotation,
                                       n_candidates)])

blood_sig <- blood_in[fdr < FDR_THRESHOLD, .(
  compartment, timepoint, feature = toupper(feature), mz, rt, ion_mode,
  effect_size, raw_p, fdr,
  direction = fifelse(effect_size < 0, "Downregulated", "Upregulated"))]
blood_sig <- merge(blood_sig, blood_ann_key,
                   by = c("compartment", "timepoint", "feature"), all.x = TRUE)

# Milk: MISAME-III BEP untargeted FDR-significant features, already annotated.
milk_ann <- fread(paste0(root, "results/milk_fdr_sig_putative_annotation.csv"))
milk_sig <- milk_ann[study == "Misame" & contrast == "BEP", .(
  compartment = "Milk",
  timepoint = fcase(visit == "14-21 days", "1421d",
                    visit == "1-2 mo.",    "pn12",
                    visit == "3-4 mo.",    "pn34",
                    default = visit),
  feature = toupper(feature), mz, rt = NA_real_, ion_mode = mode,
  effect_size = est, raw_p = NA_real_, fdr = q,
  direction = fifelse(direction == "down", "Downregulated", "Upregulated"),
  putative_name = best_annotation, n_candidates)]
milk_sig <- milk_sig[timepoint %in% c("1421d", "pn12", "pn34")]

compartment_sig <- rbindlist(list(blood_sig, milk_sig), use.names = TRUE)
# normalize the putative name: MetaboAnalyst-style "- none -"/blank -> NA
compartment_sig[putative_name %in% c("", "- none -", "none", "NA"), putative_name := NA_character_]
compartment_sig[, normalized_name := gsub("[^a-z0-9]+", "", tolower(putative_name))]
compartment_sig[normalized_name == "", normalized_name := NA_character_]

# ---------------------------------------------------------------------------
# (8) 5-level supplement_status for every FDR-significant feature
# ---------------------------------------------------------------------------
# Exact feature-id match first (features on the same rLC catalogue as the APSE
# supplement), else the nearest same-ion-mode supplement feature within 25 ppm.
compartment_sig[, row_id := .I]
prof <- supplement_profile[!duplicated(feature),
                           .(feature, supplement_feature, supplement_mz, ion_mode,
                             supplement_detected, supplement_reliably_detected,
                             supplement_abundant)]

# tier 1: exact feature id
exact <- merge(compartment_sig[, .(row_id, feature)], prof, by = "feature")[
  , .(row_id, supplement_feature, supplement_detected, supplement_reliably_detected,
      supplement_abundant, match_type = "Exact feature identifier", supplement_match_ppm = 0)]

# tier 2: nearest same-ion-mode supplement feature within 25 ppm, for the rest
mass_src <- compartment_sig[!(row_id %in% exact$row_id) & !is.na(mz) & !is.na(ion_mode),
                            .(row_id, mz = as.numeric(mz), mode = ion_mode)]
mass <- exact[0]                                     # empty, same schema
if (nrow(mass_src)) {
  A <- copy(mass_src); B <- prof[!is.na(supplement_mz),
        .(supplement_feature, mz = supplement_mz, mode = ion_mode,
          supplement_detected, supplement_reliably_detected, supplement_abundant)]
  A[, `:=`(mz_lo = mz * (1 - MASS_TOL_PPM / 1e6), mz_hi = mz * (1 + MASS_TOL_PPM / 1e6))]
  B[, `:=`(mz_lo = mz, mz_hi = mz)]
  setkey(A, mz_lo, mz_hi); setkey(B, mz_lo, mz_hi)
  h <- foverlaps(B, A, type = "within", nomatch = 0L)
  h <- h[!is.na(mode) & !is.na(i.mode) & mode == i.mode]
  if (nrow(h)) {
    # foverlaps(x=B, y=A): columns unique to B (supplement_*) stay unprefixed;
    # only the shared mz/mode clash, where A's copy takes the `i.` prefix.
    h[, ppm := abs(mz - i.mz) / i.mz * 1e6]
    setorder(h, ppm)
    h <- h[!duplicated(row_id)]                      # nearest supplement feature per sig feature
    mass <- h[, .(row_id, supplement_feature, supplement_detected,
                  supplement_reliably_detected, supplement_abundant,
                  match_type = "Mass within 25 ppm", supplement_match_ppm = ppm)]
  }
}

status_tbl <- rbindlist(list(exact, mass), use.names = TRUE)
compartment_sig <- merge(compartment_sig, status_tbl, by = "row_id", all.x = TRUE)

compartment_sig[, supplement_status := fcase(
  is.na(supplement_feature),                            "No matched supplement feature",
  supplement_abundant == TRUE,                          "Supplement-abundant",
  supplement_reliably_detected == TRUE,                 "Reliably detected",
  supplement_detected == TRUE,                          "Detected",
  default = "Not detected")]

fwrite(compartment_sig, paste0(root, "results/supplement_status_fdr_features.csv"))

status_summary <- compartment_sig[, .(n_features = .N),
                                  by = .(compartment, timepoint, direction, supplement_status)]
cat("\nsupplement status of FDR-significant features:\n")
print(status_summary[order(compartment, timepoint, direction, supplement_status)])

# ---------------------------------------------------------------------------
# (9) Putative-name cross-compartment link tier (direction-only concordance)
# ---------------------------------------------------------------------------
# Features linked across DIFFERENT compartments by an identical normalized
# putative name. Concordance is by direction only (intensities are not comparable
# across platforms). This is the name-based tier the repo's id/mass matcher lacked.
named <- compartment_sig[!is.na(normalized_name)]
name_pairs <- merge(
  named[, .(normalized_name, c1 = compartment, t1 = timepoint, f1 = feature,
            dir1 = direction, name1 = putative_name)],
  named[, .(normalized_name, c2 = compartment, t2 = timepoint, f2 = feature,
            dir2 = direction)],
  by = "normalized_name", allow.cartesian = TRUE)[as.character(c1) < as.character(c2)]
name_pairs[, same_direction := dir1 == dir2]
name_concordance <- name_pairs[, .(n_links = .N), by = .(c1, c2, same_direction)][order(c1, c2)]
cat("\nputative-name cross-compartment links (direction concordance):\n")
print(name_concordance)

# ---------------------------------------------------------------------------
# (10a) Stacked supplement-status bar
# ---------------------------------------------------------------------------
ct_levels <- c("Maternal plasma · incl", "Maternal plasma · tri3", "Maternal plasma · pn12",
               "Maternal VAMS · tri3", "Maternal VAMS · pn56",
               "Milk · 1421d", "Milk · pn12", "Milk · pn34",
               "Infant VAMS · acco", "Infant VAMS · pn12", "Infant VAMS · pn34", "Infant VAMS · pn56")
ct_labels <- c("Maternal plasma · Enrollment", "Maternal plasma · Third trimester",
               "Maternal plasma · 1–2 months", "Maternal VAMS · Third trimester",
               "Maternal VAMS · 5–6 months", "Milk · 14–21 days", "Milk · 1–2 months",
               "Milk · 3–4 months", "Infant VAMS · Birth", "Infant VAMS · 1–2 months",
               "Infant VAMS · 3–4 months", "Infant VAMS · 5–6 months")
status_levels <- c("No matched supplement feature", "Not detected", "Detected",
                   "Reliably detected", "Supplement-abundant")
status_colors <- c("No matched supplement feature" = "#D9D9D9", "Not detected" = "#969696",
                   "Detected" = "#F2CC8F", "Reliably detected" = "#59A89C",
                   "Supplement-abundant" = "#7B2CBF")

bar_df <- copy(status_summary)
bar_df[, compartment_timepoint := factor(paste(compartment, timepoint, sep = " · "),
                                         levels = ct_levels, labels = ct_labels)]
bar_df[, supplement_status := factor(supplement_status, levels = status_levels)]
bar_df[, direction := factor(direction, levels = c("Upregulated", "Downregulated"))]
bar_df <- bar_df[!is.na(compartment_timepoint)]

supplement_status_plot <- ggplot(bar_df,
    aes(x = n_features, y = compartment_timepoint, fill = supplement_status)) +
  geom_col(position = "fill", width = 0.75, color = "white", linewidth = 0.25) +
  facet_wrap(vars(direction), ncol = 1) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = status_colors, drop = FALSE) +
  labs(x = "Percentage of FDR-significant features", y = NULL, fill = "BEP supplement status") +
  theme_bw(base_size = 8) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.box = "vertical") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave(paste0(root, "figures/supplement_status_by_compartment.png"),
       supplement_status_plot, width = 12, height = 8, units = "in", dpi = 600, bg = "white")
cat("\nwrote figures/supplement_status_by_compartment.png\n")

# ---------------------------------------------------------------------------
# (10b) Supplement-status-shaded cross-compartment tracking line plot
# ---------------------------------------------------------------------------
# Track putative metabolites that recur in >= 2 compartment x timepoint positions
# (name-based linkage), one line per metabolite, points shaped by supplement status.
track_df <- compartment_sig[!is.na(putative_name)]
track_df[, compartment_timepoint := factor(paste(compartment, timepoint, sep = " · "),
                                           levels = ct_levels, labels = ct_labels)]
track_df <- track_df[!is.na(compartment_timepoint)]
track_df[, tracking_name := putative_name]
track_df[, supplement_label := fcase(
  supplement_abundant == TRUE,          "Supplement-abundant",
  supplement_reliably_detected == TRUE, "Detected in supplement",
  default = "Not reliably detected in supplement")]
# keep one row per metabolite x position (most significant), and metabolites seen
# in at least two distinct positions
setorder(track_df, direction, tracking_name, compartment_timepoint, fdr)
track_df <- track_df[!duplicated(paste(direction, tracking_name, compartment_timepoint))]
track_df[, n_positions := uniqueN(compartment_timepoint), by = .(direction, tracking_name)]
track_df <- track_df[n_positions >= 2]

if (nrow(track_df)) {
  track_labels <- track_df[, .SD[which.min(as.integer(compartment_timepoint))],
                           by = .(direction, tracking_name)]
  compartment_tracking_plot <- ggplot(track_df,
      aes(x = compartment_timepoint, y = effect_size, group = tracking_name)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_line(color = "grey70", alpha = 0.7, linewidth = 0.5) +
    geom_point(aes(shape = supplement_label, color = direction), size = 2.5, alpha = 0.9) +
    geom_text_repel(data = track_labels, aes(label = tracking_name, color = direction),
                    size = 2.3, box.padding = 0.5, min.segment.length = 0,
                    max.overlaps = Inf, seed = 123, show.legend = FALSE) +
    scale_shape_manual(values = c("Supplement-abundant" = 17, "Detected in supplement" = 15,
                                  "Not reliably detected in supplement" = 16)) +
    scale_color_manual(values = c("Upregulated" = "#E15759", "Downregulated" = "#4E79A7")) +
    labs(x = NULL, y = "Average treatment effect", shape = "BEP supplement", color = "Direction") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

  ggsave(paste0(root, "figures/compartment_tracking_supplement_status.png"),
         compartment_tracking_plot, width = 14, height = 9, units = "in", dpi = 600, bg = "white")
  cat(sprintf("wrote figures/compartment_tracking_supplement_status.png (%d tracked metabolite-positions)\n",
              nrow(track_df)))
} else {
  cat("no metabolites recur across >= 2 compartment x timepoint positions; tracking plot skipped\n")
}

# ---------------------------------------------------------------------------
# Export consolidated workbook
# ---------------------------------------------------------------------------
writexl::write_xlsx(list(
  supplement_profile   = as.data.frame(supplement_profile),
  supplement_status    = as.data.frame(compartment_sig),
  status_summary       = as.data.frame(status_summary),
  name_concordance     = as.data.frame(name_concordance),
  count_validation     = as.data.frame(count_validation)
), paste0(root, "results/compartment_tracking_supplement.xlsx"))
cat("wrote results/compartment_tracking_supplement.xlsx\n")
