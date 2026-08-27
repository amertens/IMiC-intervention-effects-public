# =============================================================================
# compartment-tracking.R
# Faithful R port of the cross-compartment LINKAGE at the heart of Trenton's
# "Compartment Tracking.Rmd" -- the logic that decides which metabolite features
# in different compartments (Maternal plasma, Maternal VAMS, Milk, Infant VAMS)
# are "the same metabolite" and therefore transfer.
#
# His rule (Compartment Tracking.Rmd, `compartment_pairs`, lines ~1531-1629), with
# `mass_tolerance_ppm <- 25`:
#   Two SIGNIFICANT features in DIFFERENT compartments are linked if ANY of:
#     (1) exact feature identifier            feature_1 == feature_2
#     (2) identical normalized putative name   normalized_name_1 == normalized_name_2
#     (3) same ionization mode AND mass within 25 ppm
#   Direction concordance (same_direction) is tracked; cross-compartment transfer
#   uses direction only (absolute intensities are not comparable across platforms).
#   normalized_name = str_to_lower(putative_name) with [^a-z0-9]+ stripped.
#
# WHY THIS PORT EXISTS: Panel D previously linked features by m/z rounded to 3
# decimals (~1.7 ppm) -- ~15x tighter than his 25 ppm window and with NO name /
# feature-id matching. That split features his method merges: e.g. the headline
# carnitine (m/z 286.198 in milk/plasma vs 286.202 in infant/VAMS, ~14 ppm apart)
# was broken into two 2-3 compartment components, so nothing reached 4 compartments
# and the milk->infant transfer disappeared. This port reproduces his linkage and
# is the source of truth Panel D now consumes (fig6D-crosscompartment.R).
#
# SCOPE: this ports the LINKAGE + concordant-transfer identification, reproducibly
# from the in-repo slim results. The upstream mummichog putative-name annotation
# and the APSE-supplement 25-ppm matching in his Rmd depend on his external
# mummichog result folders / supplement feature list (already folded into the slim
# RDS `putative_name`, and handled by the annotation hand-off); they are not
# re-run here.
#
# Run from repo root:  Rscript src/trenton-ports/compartment-tracking.R
# =============================================================================
suppressMessages({ library(dplyr); library(stringr); library(tidyr) })

SLIM <- "results/compartment_tracking/compartment_named_features_slim.RDS"
OUT  <- "results/trenton-ports/compartment_tracking"
COMP_LEVELS <- c("Maternal plasma", "Maternal VAMS", "Milk", "Infant VAMS")
MASS_TOLERANCE_PPM <- 25   # his `mass_tolerance_ppm`

# his normalized_name: lower-case, strip everything but a-z0-9
ct_normalized_name <- function(putative_name) {
  str_replace_all(str_to_lower(putative_name), "[^a-z0-9]+", "")
}

# robust significant-flag coercion (logical / "TRUE" / 1 / factor)
ct_as_sig <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x == 1)
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "significant")
}

# ---------------------------------------------------------------------------
# ct_compartment_pairs(): his exact pairwise linkage over SIGNIFICANT features.
# Input must carry: compartment, timepoint, feature, mz, ion_mode, direction,
# effect_size, fdr, putative_name. Returns matched cross-compartment pairs with a
# match_type and same_direction flag.
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
# ct_tracking_components(): union-find over SAME-DIRECTION matched pairs, so a
# "compound" is a connected component of features his rule judges to be the same
# metabolite moving in one direction. Returns the significant features with a
# `tracking_group` id (features in no same-direction cross-compartment pair get
# their own singleton group). This generalises his linked_{up,down}_ids to the
# >=2-compartment concordant-transfer set Panel D needs.
# ---------------------------------------------------------------------------
ct_tracking_components <- function(sig) {
  sig <- sig %>%
    mutate(compartment = as.character(compartment), timepoint = as.character(timepoint),
           tracking_id = paste(compartment, timepoint, feature, sep = "__"))
  pairs <- ct_compartment_pairs(sig) %>% filter(same_direction)
  parent <- setNames(sig$tracking_id, sig$tracking_id)
  find <- function(x) { while (parent[[x]] != x) { parent[[x]] <<- parent[[ parent[[x]] ]]; x <- parent[[x]] }; x }
  union <- function(x, y) { rx <- find(x); ry <- find(y); if (rx != ry) parent[[rx]] <<- ry }
  if (nrow(pairs)) for (i in seq_len(nrow(pairs))) union(pairs$tracking_id_1[i], pairs$tracking_id_2[i])
  sig$tracking_group <- vapply(sig$tracking_id, find, character(1))
  sig
}

# ---------------------------------------------------------------------------
# ct_trenton_track(): Trenton's NAME-BASED grouping (Compartment Tracking.Rmd
# 1743-1749, 2556-2562) -- tracking_name = the feature's own putative name (or
# "m/z <round4>" if none), grouped WITHIN `direction`. Unlike ct_tracking_components()
# above, this does NOT merge isobars/near-mass matches into one component -- each
# feature keeps its own annotation as its group label. Adds `tracking_name` plus
# per-group `n_cells` (distinct compartment x timepoint cells) and `n_comp`
# (distinct compartments) to every row; filtering is left to the caller:
#   - Fig S10 keeps n_cells >= 2 (his literal rule, lines 2556-2562)
#   - Fig 6D additionally requires n_comp >= 2 (genuine cross-COMPARTMENT
#     transfer -- what its "N compartment-level estimates" legend describes)
# This is the single source of the name-based rule so Fig 6D, Fig S10 and Table
# S11 cannot drift onto different groupings again -- see
# Manuscript/PROVENANCE_AUDIT.md §6.4.
# ---------------------------------------------------------------------------
ct_trenton_track <- function(sig) {
  sig %>%
    mutate(compartment = as.character(compartment), timepoint = as.character(timepoint),
           tracking_name = dplyr::coalesce(putative_name, paste0("m/z ", round(mz, 4))),
           .ct_cell = paste(compartment, timepoint, sep = " · ")) %>%
    group_by(direction, tracking_name) %>%
    mutate(n_cells = dplyr::n_distinct(.ct_cell), n_comp = dplyr::n_distinct(compartment)) %>%
    ungroup() %>% select(-.ct_cell)
}

if (sys.nframe() == 0) {
  d <- readRDS(SLIM) %>% filter(!is.na(mz), compartment %in% COMP_LEVELS)
  sig <- d %>% filter(ct_as_sig(significant))
  message("significant compartment features: ", nrow(sig))

  pairs <- ct_compartment_pairs(sig)
  comp  <- ct_tracking_components(sig)

  # concordant transfer: a tracking_group significant in >= 2 distinct compartments
  transfer <- comp %>% group_by(tracking_group) %>%
    summarise(n_comp = n_distinct(compartment), direction = first(direction),
              compartments = paste(sort(unique(compartment)), collapse = "+"),
              name = { x <- putative_name[!is.na(putative_name) & putative_name != ""]; if (length(x)) x[1] else NA_character_ },
              .groups = "drop") %>%
    filter(n_comp >= 2) %>% arrange(desc(n_comp))

  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  write.csv(pairs,    file.path(OUT, "compartment_pairs.csv"), row.names = FALSE, na = "")
  write.csv(transfer, file.path(OUT, "concordant_transfer_components.csv"), row.names = FALSE, na = "")
  message(sprintf("linked pairs: %d | concordant-transfer components (>=2 compartments): %d | reaching 4: %d",
                  nrow(pairs), nrow(transfer), sum(transfer$n_comp == 4)))
  message("-> ", OUT)
}
