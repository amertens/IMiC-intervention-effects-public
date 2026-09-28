# =============================================================================
# compartment-tracking.R
#
# Helper functions for the cross-compartment linkage behind Fig 6D and Table S11,
# ported from the original R Markdown compartment-tracking analysis. They decide
# which significant metabolite features in different compartments (maternal
# plasma, maternal VAMS, milk, infant VAMS) count as the same metabolite, and
# group features by name. Sourced by src/metaboanalyst/build-compartment-query-list.R
# and src/metaboanalyst/run-compartment-pathway.R; defines functions only.
#
# Linkage rule (mass tolerance 25 ppm): two significant features in different
# compartments are linked if any of
#   (1) the feature identifiers are equal,
#   (2) the normalized putative names are equal (lower case, non-alphanumerics removed), or
#   (3) the ion modes are equal and the masses are within 25 ppm.
# Transfer is judged on direction only, because absolute intensities are not
# comparable across platforms. The putative names themselves come from the
# upstream mummichog annotation already present in the input table; that
# annotation is not re-run here.
#
# Inputs : none (functions take a data frame of significant features)
# Outputs: none
# =============================================================================
suppressMessages({ library(dplyr); library(stringr); library(tidyr) })

MASS_TOLERANCE_PPM <- 25

# lower-case, strip everything but a-z0-9
ct_normalized_name <- function(putative_name) {
  str_replace_all(str_to_lower(putative_name), "[^a-z0-9]+", "")
}

# ---------------------------------------------------------------------------
# ct_compartment_pairs(): pairwise linkage over significant features.
# Input must carry: compartment, timepoint, feature, mz, ion_mode, direction,
# effect_size, fdr, putative_name. Returns matched cross-compartment pairs with a
# match_type and same_direction flag; each pair appears once, keeping the best
# match (exact identifier first, then smallest ppm difference, then lowest FDR).
# ---------------------------------------------------------------------------
ct_compartment_pairs <- function(sig) {
  sig <- sig %>%
    mutate(compartment = as.character(compartment), timepoint = as.character(timepoint),
           normalized_name = ct_normalized_name(putative_name),
           tracking_id = paste(compartment, timepoint, feature, sep = "__"))
  keep <- c("tracking_id","compartment","timepoint","ion_mode","feature","mz",
            "effect_size","direction","fdr","putative_name","normalized_name")
  a <- sig %>% select(all_of(keep)) %>% rename_with(~ paste0(.x, "_1"))
  b <- sig %>% select(all_of(keep)) %>% rename_with(~ paste0(.x, "_2"))
  tidyr::crossing(a, b) %>%
    filter(compartment_1 < compartment_2) %>%                       # cross-compartment, ordered
    mutate(
      ppm_difference      = abs(mz_1 - mz_2) / ((mz_1 + mz_2) / 2) * 1e6,
      exact_feature_match = feature_1 == feature_2,
      mass_match          = ion_mode_1 == ion_mode_2 & ppm_difference <= MASS_TOLERANCE_PPM,
      similar_name_match  = !is.na(normalized_name_1) & !is.na(normalized_name_2) &
                            normalized_name_1 != "" & normalized_name_1 == normalized_name_2,
      same_direction      = direction_1 == direction_2,
      match_type = case_when(exact_feature_match ~ "Exact feature identifier",
                             similar_name_match  ~ "Same putative name",
                             mass_match          ~ "Mass within 25 ppm",
                             TRUE ~ NA_character_)) %>%
    filter(!is.na(match_type)) %>%
    arrange(desc(exact_feature_match), ppm_difference, fdr_1, fdr_2) %>%
    distinct(tracking_id_1, tracking_id_2, .keep_all = TRUE)
}

# ---------------------------------------------------------------------------
# ct_name_track(): name-based grouping. tracking_name = the feature's own
# putative name (or "m/z <mz rounded to 4 dp>" if none), grouped within
# `direction`. Isobars and near-mass matches are not merged: each feature keeps
# its own annotation as its group label. Adds `tracking_name` plus per-group
# `n_cells` (distinct compartment x timepoint cells) and `n_comp` (distinct
# compartments) to every row; the caller chooses the threshold:
#   - the Table S11 query list uses tracking_name only and requires >= 2 cells
#     among the linked up-regulated features (build-compartment-query-list.R)
#   - Fig 6D keeps n_comp >= 2, i.e. transfer across compartments
#     (run-compartment-pathway.R)
# Both exhibits use this one function so they share the same grouping.
# ---------------------------------------------------------------------------
ct_name_track <- function(sig) {
  sig %>%
    mutate(compartment = as.character(compartment), timepoint = as.character(timepoint),
           tracking_name = dplyr::coalesce(putative_name, paste0("m/z ", round(mz, 4))),
           .ct_cell = paste(compartment, timepoint, sep = " · ")) %>%
    group_by(direction, tracking_name) %>%
    mutate(n_cells = dplyr::n_distinct(.ct_cell), n_comp = dplyr::n_distinct(compartment)) %>%
    ungroup() %>% select(-.ct_cell)
}
