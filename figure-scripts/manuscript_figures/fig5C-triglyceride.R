# =============================================================================
# fig5C-triglyceride.R
#
# Builds Fig 5C, the triglyceride -> fatty-acid composition volcano, and writes the
# fatty-acid composition table behind Table S4. The PNG is embedded as Panel C by
# fig5-tertiary-composite.R. The numerics follow the original R Markdown
# triglyceride analysis (its analyze_triglyceride_comparison step):
#   1. Read the intervention-effect results (Triglycerides, ATE only) -> sigTgs.
#      Default = combined-arm framing (one contrast per study x visit), as in the
#      figure. Set IMIC_TG_ARM_FRAMING=stratified for the per-arm contrasts.
#   2. Read the Biocrates Quant500 structure file to map each triglyceride
#      Short Name -> its LSSN molecular formula, split the LSSN into the three
#      fatty-acyl chains (sn1/sn2/sn3, with the exact string cleanups),
#      and compute each fatty acid's PROPORTION within every TG.
#   3. For every study x collection-time x contrast x DIRECTION cell: within that
#      direction (up = est>0, down = est<0), per fatty acid run a Wilcoxon test of
#      its proportion in FDR-significant vs non-significant (DIRECTION-MATCHED)
#      TGs. Both directions are analysed as explicit, separate comparisons.
#      Accumulate into trig_all.
#   4. Study-coloured volcano of the SIGNED enrichment ratio (positive = enriched
#      among UP-regulated TGs, negative = enriched among DOWN-regulated TGs) vs
#      -log10(p_value), matching Panel B's conventions.
#
# Inputs:  results/combined_intervention_effects_results_combined_arms.RDS
#            (or ..._stratified_arms.RDS; src/2 analysis/clean_results.R)
#          data/additional datasets/Copy of Quant 500_BioIDs_20191031.xlsx, the
#            Biocrates Quant 500 "BioIDs" structure file (sheet "FIA Part"), or the
#            path in IMIC_BIOCRATES_BIOIDS
# Outputs: figures/figure5_panelC_tg_composition{,_nolegend}.png
#            (_stratified{,_nolegend} for the stratified framing)
#          results/metaboanalyst/triglyceride_fa/triglyceride_fa_composition_{combined,stratified}.csv
#            (Table S4)
# [needs an on-request file] the combined-results RDS is not shipped (size).
# The BioIDs structure file is a Biocrates vendor annotation file that cannot be
# redistributed; without it the triglycerides cannot be mapped to their fatty-acid
# chains, so the script stops with a message naming the expected path.
# Run from the repo root:
#   Rscript "figure-scripts/manuscript_figures/fig5C-triglyceride.R"
# =============================================================================

suppressMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
})
source("figure-scripts/manuscript_figures/study_colors.R")   # shared study colours and shapes
source("figure-scripts/0_figure-functions.R")                # theme_imic(), imic_logp_title

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# The figure uses the combined-arm framing; "stratified" gives the per-arm contrasts
# used for the MISAME-III rows of Table S4. Set via env var so it can be re-run
# without editing this file: IMIC_TG_ARM_FRAMING=stratified Rscript ...fig5C-triglyceride.R
ARM_FRAMING <- Sys.getenv("IMIC_TG_ARM_FRAMING", "combined")   # "combined" | "stratified"

RESULTS_RDS <- if (identical(ARM_FRAMING, "stratified")) {
  "results/combined_intervention_effects_results_stratified_arms.RDS"
} else {
  "results/combined_intervention_effects_results_combined_arms.RDS"
}

# Biocrates Quant 500 structure file. Expected location plus an env-var override
# for a copy elsewhere.
BIOCRATES_BIOIDS_PATH <- Sys.getenv(
  "IMIC_BIOCRATES_BIOIDS",
  "data/additional datasets/Copy of Quant 500_BioIDs_20191031.xlsx"
)

# The printed panel is the combined-arm framing; a stratified run (for the Table S4
# MISAME-III rows) must not overwrite it, so it writes to a distinctly named PNG.
OUT_PNG <- if (identical(ARM_FRAMING, "stratified")) {
  "figures/figure5_panelC_tg_composition_stratified.png"
} else {
  "figures/figure5_panelC_tg_composition.png"
}
# No-legend twin, embedded by the Fig 5 composite (fig5-tertiary-composite.R),
# which draws one shared Study legend under Panels B+C instead of a per-panel one.
OUT_PNG_NOLEGEND <- if (identical(ARM_FRAMING, "stratified")) {
  "figures/figure5_panelC_tg_composition_stratified_nolegend.png"
} else {
  "figures/figure5_panelC_tg_composition_nolegend.png"
}

# ---------------------------------------------------------------------------
# Stop with a clear message if the Biocrates structure file (a vendor file that
# is not in the repository) is missing.
# ---------------------------------------------------------------------------
.biocrates_reference <- function(path = BIOCRATES_BIOIDS_PATH) {
  if (is.null(path) || !file.exists(path)) {
    stop(
      "Biocrates Quant500 BioIDs structure file not found at '", path, "'.\n",
      "This file ('Copy of Quant 500_BioIDs_20191031.xlsx', sheet 'FIA Part') ",
      "is required to map each triglyceride to its constituent fatty-acid ",
      "chains; without it the fatty-acid composition volcano cannot be built. ",
      "It is a Biocrates vendor annotation file that cannot be redistributed, ",
      "so it is not stored in the repository. Place it at the path above, ",
      "or set the IMIC_BIOCRATES_BIOIDS environment variable to its location.",
      call. = FALSE
    )
  }
  path
}

# ---------------------------------------------------------------------------
# Step 1. Intervention-effect results -> significant triglycerides (ATE).
# ---------------------------------------------------------------------------
load_sig_tgs <- function(rds_path = RESULTS_RDS) {
  readRDS(rds_path) %>%
    filter(category == "Triglycerides", measure == "ATE")
}

# ---------------------------------------------------------------------------
# Step 2. Biocrates structure -> fatty-acid proportions within each TG.
#   As in the original R Markdown analysis (all sn1/sn2/sn3 str_replace cleanups kept).
# ---------------------------------------------------------------------------
build_proportional_data <- function(bioids_path) {
  # Read the "FIA Part" sheet (guard already checked the file exists).
  rawData <- read_excel(bioids_path, sheet = "FIA Part")

  data <-
    rawData %>%
    rename(shortName = `Short Name`, lssn = `Molecule`) %>%
    filter(str_detect(shortName, "^TG")) %>%
    separate(lssn, into = c("sn1", "sn2", "sn3"), sep = "/",
             remove = FALSE, fill = "right") %>%
    mutate(
      sn1 = str_replace(sn1, "^TG\\(", ""),
      sn1 = str_replace(sn1, "(\\d)\\)$", "\\1"),
      sn1 = if_else(str_starts(sn1, "PA") | str_starts(sn1, "PG"), paste0(sn1, ")"), sn1),
      sn1 = replace_na(sn1, "unknown"),
      sn2 = str_replace(sn2, "^TG\\(", ""),
      sn2 = str_replace(sn2, "(\\d)\\)$", "\\1"),
      sn2 = if_else(str_starts(sn2, "PA") | str_starts(sn2, "PG"), paste0(sn2, ")"), sn2),
      sn2 = replace_na(sn2, "unknown"),
      sn3 = str_replace(sn3, "^TG\\(", ""),
      sn3 = str_replace(sn3, "(\\d)\\)$", "\\1"),
      sn3 = if_else(str_starts(sn3, "PA") | str_starts(sn3, "PG"), paste0(sn3, ")"), sn3),
      sn3 = replace_na(sn3, "unknown"),
      sn3 = str_remove(sn3, "\\)\\[iso\\d+\\]"),
      sn3 = if_else(str_detect(sn3, "\\)\\)$"), str_replace(sn3, "\\)\\)$", ")"), sn3),
      sn3 = if_else(str_detect(sn3, ":0\\)$"), str_replace(sn3, "\\)$", ""), sn3),
      sn1 = if_else(str_detect(sn1, "19:1\\(9Z\\)\\)$"), str_replace(sn1, "19:1\\(9Z\\)\\)$", "19:1(9Z)"), sn1),
      sn1 = if_else(str_detect(sn1, "20:1\\(11Z\\)\\)$"), str_replace(sn1, "20:1\\(11Z\\)\\)$", "20:1(11Z)"), sn1),
      sn1 = if_else(str_detect(sn1, "19:0\\)$"), str_replace(sn1, "19:0\\)$", "19:0"), sn1),
      sn1 = str_replace(sn1, "PA\\(19:1\\(9Z\\)\\)$", "PA(19:1(9Z)"),
      sn1 = str_replace(sn1, "PA\\(20:1\\(11Z\\)\\)$", "PA(20:1(11Z)"),
      sn2 = if_else(str_detect(sn2, "19:1\\(9Z\\)\\)$"), str_replace(sn2, "19:1\\(9Z\\)\\)$", "19:1(9Z)"), sn2),
      sn2 = if_else(str_detect(sn2, "20:1\\(11Z\\)\\)$"), str_replace(sn2, "20:1\\(11Z\\)\\)$", "20:1(11Z)"), sn2),
      sn2 = if_else(str_detect(sn2, "19:0\\)$"), str_replace(sn2, "19:0\\)$", "19:0"), sn2),
      sn3 = if_else(str_detect(sn3, "19:1\\(9Z\\)\\)$"), str_replace(sn3, "19:1\\(9Z\\)\\)$", "19:1(9Z)"), sn3),
      sn3 = if_else(str_detect(sn3, "20:1\\(11Z\\)\\)$"), str_replace(sn3, "20:1\\(11Z\\)\\)$", "20:1(11Z)"), sn3),
      sn3 = if_else(str_detect(sn3, "19:0\\)$"), str_replace(sn3, "19:0\\)$", "19:0"), sn3),
      sn3 = if_else(str_detect(sn3, "20:2n6\\)$"), str_replace(sn3, "20:2n6\\)$", "20:2n6"), sn3),
      sn3 = if_else(str_detect(sn3, "20:3n6\\)$"), str_replace(sn3, "20:3n6\\)$", "20:3n6"), sn3)
    ) %>%
    select(shortName, lssn, sn1, sn2, sn3)

  # Drop TGs whose shortName == LSSN (unresolved composition), pivot the three
  # acyl positions long, and compute each fatty acid's proportion within its TG.
  longfilteredData <-
    data %>%
    filter(shortName != lssn) %>%
    pivot_longer(cols = c(sn1, sn2, sn3),
                 names_to = "snPosition", values_to = "fattyAcid") %>%
    select(-c(snPosition, lssn))

  proportionalData <-
    longfilteredData %>%
    group_by(shortName, fattyAcid) %>%
    tally() %>%
    group_by(shortName) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup() %>%
    select(-n)

  proportionalData
}

# ---------------------------------------------------------------------------
# Step 3. Per-cell fatty-acid composition test.
#
# One helper applied over the study x collection-time x contrast x DIRECTION grid:
#   * Both directions are analysed as explicit, separate comparisons.
#   * The non-significant REFERENCE group is DIRECTION-MATCHED: significant vs
#     non-significant TGs are compared only within the same effect direction (all
#     TGs in the cell already share est>0 for "up" / est<0 for "down").
#   * mean_diff = mean(proportion | FDR-significant) - mean(proportion | non-sig);
#     p_value from wilcox.test(proportion ~ sigFDR) per fatty acid, only where BOTH
#     groups are present.
# Numeric definitions match the original analysis (analyze_triglyceride_comparison).
# ---------------------------------------------------------------------------
analyze_cell <- function(studytime_val, contrast_val, direction_val,
                         sigTgs, proportionalData) {

  tg_list <-
    sigTgs %>%
    filter(studytime == studytime_val, contrast == contrast_val,
           if (direction_val == "up") est > 0 else est < 0) %>%
    mutate(shortName = str_replace(label, "Triacylglyceride", "TG"),
           shortName = str_replace(shortName, "TG \\(", "TG(")) %>%
    ungroup() %>%
    select(shortName, est, sigFDR, pval_adj) %>%
    filter(!is.na(pval_adj), sigFDR %in% c(0, 1)) %>%
    distinct() %>%
    inner_join(proportionalData, by = "shortName")

  # A cell contributes nothing unless it has both a non-significant (sigFDR==0)
  # background and an FDR-significant (sigFDR==1) group WITHIN this direction;
  # this is why many cells are simply absent from the pooled result (e.g. no
  # FDR-significant down-regulated TGs in that stratum).
  if (nrow(tg_list) == 0 ||
      !any(tg_list$sigFDR == 0, na.rm = TRUE) ||
      !any(tg_list$sigFDR == 1, na.rm = TRUE)) {
    return(NULL)
  }

  # Wilcoxon per fatty acid, only where both groups are present.
  wilcox_results <-
    tg_list %>%
    group_by(fattyAcid) %>%
    filter(sum(sigFDR == 0, na.rm = TRUE) > 0,
           sum(sigFDR == 1, na.rm = TRUE) > 0) %>%
    summarise(
      test = list(wilcox.test(proportion ~ sigFDR, na.action = na.omit, exact = FALSE)),
      .groups = "drop"
    ) %>%
    mutate(p_value = map_dbl(test, "p.value")) %>%
    select(-test)

  # Mean difference in proportion (significant minus background) per fatty acid.
  mean_diff_tbl <-
    tg_list %>%
    group_by(fattyAcid, sigFDR) %>%
    summarise(mean_proportion = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sigFDR, values_from = mean_proportion,
                names_prefix = "sigFDR_")

  if (!all(c("sigFDR_0", "sigFDR_1") %in% names(mean_diff_tbl))) return(NULL)

  mean_diff_tbl %>%
    filter(!is.na(sigFDR_0), !is.na(sigFDR_1)) %>%
    mutate(mean_diff = sigFDR_1 - sigFDR_0) %>%
    inner_join(wilcox_results, by = "fattyAcid") %>%
    mutate(
      studytime = studytime_val,
      contrast  = contrast_val,
      direction = direction_val,
      # Split "Misame (14-21 days)" -> study + timepoint.
      study     = str_trim(str_replace(studytime_val, "\\s*\\(.*\\)$", "")),
      timepoint = str_replace(str_extract(studytime_val, "\\(.*\\)$"), "^\\((.*)\\)$", "\\1")
    ) %>%
    arrange(p_value)
}

# ---------------------------------------------------------------------------
# Full-name / abbreviation lookup for fatty-acid labels (from the original analysis).
# ---------------------------------------------------------------------------
annotate_fatty_acids <- function(trig_all) {
  trig_all %>%
    mutate(
      full_name_fatty_acid = case_when(
        fattyAcid == "13:0" ~ "Tridecylic Acid",
        fattyAcid == "14:0" ~ "Myristic Acid",
        fattyAcid == "14:1(9Z)" ~ "Myristoleic Acid",
        fattyAcid == "15:0" ~ "Pentadecylic Acid",
        fattyAcid == "15:1(9Z)" ~ "9Z-Pentadecenoic Acid",
        fattyAcid == "16:0" ~ "Palmitic Acid",
        fattyAcid == "16:1(9Z)" ~ "Palmitoleic Acid",
        fattyAcid == "17:0" ~ "Margaric Acid",
        fattyAcid == "19:0" ~ "Nonadecylic Acid",
        fattyAcid == "19:1(9Z)" ~ "9Z-Nonadecenoic Acid",   # 19:1 monoene (was "Nonadecanoic", the 19:0 saturated name)
        fattyAcid == "21:0" ~ "Heneicosylic Acid",
        fattyAcid == "22:2(13Z,16Z)" ~ "13Z,16Z-Docosadienoic Acid",
        fattyAcid == "22:4(7Z,10Z,13Z,16Z)" ~ "Adrenic Acid",
        fattyAcid == "22:5(7Z,10Z,13Z,16Z,19Z)" ~ "Docosapentaenoic Acid",
        fattyAcid == "18:3(9Z,12Z,15Z)" ~ "α-Linolenic Acid",   # 18:3n-3 (Δ9,12,15)
        fattyAcid == "18:3(6Z,9Z,12Z)" ~ "γ-Linolenic Acid",    # 18:3n-6 (Δ6,9,12)
        fattyAcid == "18:4(6Z,9Z,12Z,15Z)" ~ "Stearidonic Acid",
        fattyAcid == "20:2(11Z,14Z)" ~ "Eicosadienoic Acid",
        fattyAcid == "20:3(8Z,11Z,14Z)" ~ "Dihomo-γ-Linolenic Acid",
        fattyAcid == "20:3(5Z,8Z,11Z)" ~ "Mead Acid",
        fattyAcid == "20:4(5Z,8Z,11Z,14Z)" ~ "Arachidonic Acid",
        fattyAcid == "20:1(11Z)" ~ "Gondoic Acid",
        fattyAcid == "12:0" ~ "Lauric Acid",
        fattyAcid == "17:2(9Z,12Z)" ~ "9,12-Heptadecadienoic Acid",   # 17:2 diene (was "Margaric acid", the 17:0 saturated name)
        fattyAcid == "20:2n6" ~ "11,14-Eicosadienoic Acid",           # full name for 20:2n6 (previously NA; abbr 11,14-EDA only)
        TRUE ~ NA_character_
      ),
      abbr = case_when(
        fattyAcid == "18:3(9Z,12Z,15Z)" ~ "ALA",   # 18:3n-3
        fattyAcid == "18:3(6Z,9Z,12Z)" ~ "GLA",    # 18:3n-6
        fattyAcid == "18:4(6Z,9Z,12Z,15Z)" ~ "SDA",
        fattyAcid == "20:2(11Z,14Z)" ~ "EDA",
        fattyAcid == "20:3(8Z,11Z,14Z)" ~ "DGLA",
        fattyAcid == "20:4(5Z,8Z,11Z,14Z)" ~ "ARA",
        fattyAcid == "22:4(7Z,10Z,13Z,16Z)" ~ "AdA",
        fattyAcid == "22:5(7Z,10Z,13Z,16Z,19Z)" ~ "DPA",
        fattyAcid == "20:2n6" ~ "11,14-EDA",
        TRUE ~ NA_character_
      )
    )
}

# ---------------------------------------------------------------------------
# Step 4. The fatty-acid composition volcano ("volcano2" in the original analysis).
# ---------------------------------------------------------------------------
build_volcano <- function(trig_all, show_legend = TRUE) {
  # Fig 5C: a study-coloured fatty-acid volcano.
  #   x = SIGNED enrichment ratio. mean_diff is the proportion difference
  #       (FDR-significant minus non-significant) WITHIN a direction; we sign it so
  #       positive = fatty acids enriched among UP-regulated TGs and negative =
  #       enriched among DOWN-regulated TGs (down comparisons negated). This mirrors
  #       Panel B's "signed enrichment ratio" convention and the Fig 5C caption.
  #   y = -log10(P). ALL fatty acids are drawn (grey non-significant cloud); those
  #       nominally significant (raw P<0.05) are coloured by study and labelled with
  #       the fatty-acid abbreviation (ALA/GLA/...) or name + timepoint. Grey P<0.05
  #       line and green Q<0.05 (BH-FDR within study x timepoint x contrast x
  #       direction, matching Panel B) line; black vertical line at 0.
  study_cols <- c(imic_study_cols, "Not significant" = "grey75")
  df <- trig_all %>%
    filter(!is.na(p_value), !is.na(mean_diff)) %>%
    group_by(study, timepoint, contrast, direction) %>%
    mutate(fdr = stats::p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(enrichment_ratio = ifelse(direction == "up", mean_diff, -mean_diff),
           study_label = dplyr::recode(study, "Elicit" = "ELICIT", "Misame" = "MISAME-III",
                                       "Vital" = "Mumta-LW", .default = study),
           lab_name    = dplyr::coalesce(abbr, full_name_fatty_acid, fattyAcid),
           is_sig      = p_value < 0.05,
           is_fdr_sig  = fdr < 0.05,
           color_group = ifelse(is_sig, study_label, "Not significant"),
           lab         = ifelse(is_sig, paste0(lab_name, " (", timepoint, ")"), NA_character_))
  q_thr  <- suppressWarnings(max(df$p_value[df$fdr < 0.05], na.rm = TRUE))
  y_q    <- if (is.finite(q_thr)) -log10(q_thr) else NA_real_
  # Only label the fatty acids discussed in the Results (allow-list on the Biocrates
  # fatty-acid code). Every point is still plotted and coloured; only the repel TEXT
  # LABELS are restricted, dropping the ~20 other nominally-significant labels
  # (Palmitoleic, Palmitic, Adrenic, Gondoic, Nonadecanoic, Margaric, Lauric, ...).
  # Labels additionally require FDR significance (is_fdr_sig), not just nominal
  # (is_sig) -- a labelled point sitting below the green Q<0.05 line was misleading.
  fa_label_allow <- c(
    "18:3(9Z,12Z,15Z)",    # alpha-Linolenic Acid (ALA)
    "18:3(6Z,9Z,12Z)",     # gamma-Linolenic Acid (GLA)
    "18:4(6Z,9Z,12Z,15Z)", # Stearidonic Acid (SDA)
    "14:0",                # Myristic Acid
    "14:1(9Z)")            # Myristoleic Acid (borderline; kept)
  lab_df <- dplyr::filter(df, is_fdr_sig, fattyAcid %in% fa_label_allow)

  ggplot(df, aes(x = enrichment_ratio, y = -log10(p_value))) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    # line colours/weights and captions as in Figs 3B and 5B
    geom_hline(yintercept = -log10(0.05), color = "#BAB0AC", linewidth = 0.4) +
    { if (is.finite(y_q)) geom_hline(yintercept = y_q, color = "#59A14F", linewidth = 0.4) } +
    # captions at the right edge (the left holds the Myristic Acid labels), as in 5B
    { if (is.finite(y_q)) annotate("text", x = Inf, y = y_q, parse = TRUE,
                                   label = 'italic(Q)*"-value" < 0.05',
                                   hjust = 1.05, vjust = -0.4, size = 2.4, color = "#3C8C3C") } +
    annotate("text", x = Inf, y = -log10(0.05), parse = TRUE, label = 'italic(P)*"-value" < 0.05',
             hjust = 1.05, vjust = -0.4, size = 2.4, color = "#8A8580") +
    geom_point(aes(color = color_group, shape = color_group), size = 2, alpha = 0.9) +  # study symbol = CVD cue
    # Boxed labels (white fill + coloured border) matching Panel B's geom_label_repel.
    geom_label_repel(data = lab_df, aes(label = lab, color = color_group), size = 2.5,
                     label.padding = 0.12, box.padding = 0.3,
                     max.overlaps = Inf, min.segment.length = 0,
                     show.legend = FALSE,
                     seed = 123) +   # draw-time placement; required for byte-reproducible Panel C
    scale_color_manual(values = study_cols, name = "Study",
                       breaks = c("ELICIT", "MISAME-III", "Mumta-LW", "Not significant")) +
    scale_shape_manual(values = imic_shapes_for(names(study_cols)), name = "Study",
                       breaks = c("ELICIT", "MISAME-III", "Mumta-LW", "Not significant")) +
    guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
    labs(x = "Enrichment Ratio", y = imic_logp_title) +
    theme_imic(base_size = 9) +   # shared theme (Helvetica, font-size floors)
    theme(legend.position = if (isTRUE(show_legend)) "bottom" else "none",
          plot.title = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3)) +
    xlim(-0.5, 0.5)
}

# ---------------------------------------------------------------------------
# Table S4 export.
#
# Writes trig_all (this script's own computed columns) to a CSV named by
# ARM_FRAMING so both variants can coexist. A pure data dump of the computed
# table; it does not affect the test or the figure. Atomic (temp file + rename)
# so a partial write never lands.
#
# The printed Table S4 combines ELICIT and Mumta-LW rows from the combined-arm run
# with MISAME-III rows from the stratified run: run once with each framing and
# union the two CSVs by study.
# ---------------------------------------------------------------------------
TRIG_FA_COMBINED_CSV <- paste0(
  "results/metaboanalyst/triglyceride_fa/triglyceride_fa_composition_",
  ARM_FRAMING, ".csv")

export_trig_all <- function(trig_all, path = TRIG_FA_COMBINED_CSV) {
  # Both framings share one column schema: enrichment_ratio is the raw
  # (unsigned) mean_diff and the `direction` column carries up/down, as
  # Table S4 pairs an "Enrichment Ratio" with a separate "Regulation" column.
  # (The figure's negative-signing of down comparisons is a plot-axis transform
  # only; the tabulated enrichment ratio stays unsigned.)
  out <- data.frame(
    study            = trig_all$study,
    timepoint        = trig_all$timepoint,
    contrast         = trig_all$contrast,
    direction        = trig_all$direction,
    fatty_acid       = trig_all$fattyAcid,
    fatty_acid_name  = trig_all$full_name_fatty_acid,
    enrichment_ratio = trig_all$mean_diff,
    raw_p            = trig_all$p_value,
    stringsAsFactors = FALSE
  )
  # BH-adjusted p per analysis unit = study x timepoint x contrast x direction
  # (each comparison independently adjusted, as in the original analysis).
  out$fdr_native <- stats::ave(out$raw_p, out$study, out$timepoint, out$contrast,
                               out$direction,
                               FUN = function(p) p.adjust(p, method = "BH"))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp", Sys.getpid())
  utils::write.csv(out, tmp, row.names = FALSE, fileEncoding = "UTF-8")
  if (file.exists(path)) file.remove(path)
  file.rename(tmp, path)
  cat("wrote", path, "\n")
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
run_panel_c <- function(write = TRUE) {
  # Stop first if the Biocrates structure file is missing, so that message is
  # what the user sees.
  bioids_path <- .biocrates_reference()

  sigTgs           <- load_sig_tgs()
  proportionalData <- build_proportional_data(bioids_path)

  # Grid = every study x collection-time x contrast present in the results,
  # crossed with both effect directions (up = est>0, down = est<0).
  grid <- sigTgs %>%
    distinct(studytime, contrast) %>%
    tidyr::crossing(direction = c("up", "down"))

  trig_all <-
    purrr::pmap_dfr(
      list(grid$studytime, grid$contrast, grid$direction),
      function(st, ct, dr)
        analyze_cell(st, ct, dr, sigTgs = sigTgs, proportionalData = proportionalData)
    ) %>%
    arrange(is.na(p_value), p_value) %>%
    annotate_fatty_acids()

  # Fatty-acid composition table for Table S4 (named by ARM_FRAMING).
  export_trig_all(trig_all)

  volcano <- build_volcano(trig_all)

  if (write) {
    # Panel size 105 x 99 mm (210/2 x 297/3), as submitted, so B and C read as
    # roughly square in the Fig 5 bottom row.
    ggsave(
      filename = OUT_PNG,
      plot     = volcano,
      device   = ragg::agg_png,
      width    = 210/2,
      height   = 297/3,
      units    = "mm",
      dpi      = 300
    )
    cat("wrote", OUT_PNG, "\n")

    ggsave(
      filename = OUT_PNG_NOLEGEND,
      plot     = build_volcano(trig_all, show_legend = FALSE),
      device   = ragg::agg_png,
      width    = 210/2,
      height   = 297/3,
      units    = "mm",
      dpi      = 300
    )
    cat("wrote", OUT_PNG_NOLEGEND, "\n")
  }

  invisible(list(trig_all = trig_all, volcano = volcano))
}

if (sys.nframe() == 0) invisible(run_panel_c(write = TRUE))
