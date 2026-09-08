# build-cells.R, turn an intervention-effects tibble into analysis "cells".
#
# A "cell" is one unit of enrichment analysis: a filtered list of query
# compounds plus the reference metabolome (background) for its outcome group.
# We split each study-time x contrast combination into two directional cells so
# that up- and down-regulated compounds are enriched separately:
#   - "down": compounds whose intervention effect estimate is negative (est < 0)
#   - "up":   compounds whose intervention effect estimate is positive (est > 0)
# Downstream (run-outcome-group.R) feeds each cell's query + reference into
# MetaboAnalystR.
suppressMessages({ library(dplyr); library(stringr) })
source("src/metaboanalyst/R/label-map.R")

# The two directions we split each study/contrast into, paired with the sign of
# the effect estimate that defines them. Kept in this order (down, then up) so
# the resulting cell list (and therefore every downstream output file) is
# emitted in a stable, reproducible order.
DIRECTION_SIGNS <- c(down = -1L, up = 1L)

# Reference metabolome for one outcome group: every distinct compound measured
# in the group EXCEPT "Total *" roll-up features (e.g. "Total vitamin B1"),
# which are composite quantities, not individual metabolites the pathway library
# can map. Returned sorted so the background list is deterministic.
build_reference_names <- function(labelled_dat) {
  labelled_dat %>%
    filter(!str_starts(metaboanalyst_label, "Total")) %>%
    distinct(metaboanalyst_label) %>%
    arrange(metaboanalyst_label) %>%
    pull(metaboanalyst_label)
}

# The query compounds for a single (study-time, contrast, direction) cell:
# significant features (pval < cutoff) whose effect estimate has the requested
# sign, again excluding "Total *" roll-ups. Returns a de-duplicated character
# vector (possibly length 0, meaning "no cell here").
build_query_names <- function(labelled_dat, studytime_i, contrast_i,
                              effect_sign, pval_cutoff) {
  significant <- labelled_dat %>%
    filter(studytime == studytime_i,
           contrast == contrast_i,
           pval < pval_cutoff,
           !str_starts(metaboanalyst_label, "Total"))

  # Keep only the requested effect direction: negative estimates for "down",
  # positive for "up".
  directional <- if (effect_sign < 0) {
    filter(significant, est < 0)
  } else {
    filter(significant, est > 0)
  }

  directional %>%
    distinct(metaboanalyst_label) %>%
    pull(metaboanalyst_label)
}

# Pathway-analysis cells (Trenton's Fig 3B "Primary Outcomes (Pathway Analysis)"
# construction): ONE non-directional cell per study-time x contrast, query =
# FDR-significant features (sigFDR == 1) with positive AND negative estimates
# analysed together (excluding "Total *" roll-ups). This deliberately differs
# from the directional pval<0.05 ORA cells below: pathway impact is a topology-
# based, direction-agnostic analysis ("which pathways are perturbed"), not an
# up/down over-representation test. Matches the query in
# "Primary Outcomes (Pathway Analysis).Rmd" (filter on studytime + sigFDR == 1).
build_cells_combined_sigfdr <- function(dat, outcome_group, measure = "ATE") {
  target_group   <- outcome_group
  target_measure <- measure

  labelled_dat <- dat %>%
    filter(outcome_group == target_group, measure == target_measure) %>%
    # Exclude "Total *" roll-up features (composite quantities, not individual
    # metabolites) on the RAW label BEFORE the synonym map runs. Doing it here is
    # required because the map renames "total vitamin B1 (expressed as thiamin)"
    # -> "Thiamine", which would otherwise slip past the later
    # `!str_starts(metaboanalyst_label, "Total")` guard (and that guard is also
    # capitalization-fragile). Case-insensitive so "Total"/"total" behave alike.
    filter(!str_detect(label_f, regex("^\\s*total\\b", ignore_case = TRUE))) %>%
    mutate(metaboanalyst_label = apply_label_map(label_f))

  reference_names <- build_reference_names(labelled_dat)

  study_contrast_grid <- labelled_dat %>%
    distinct(studytime, contrast) %>%
    arrange(studytime, contrast)

  cells <- list()
  for (i in seq_len(nrow(study_contrast_grid))) {
    studytime_i <- study_contrast_grid$studytime[i]
    contrast_i  <- study_contrast_grid$contrast[i]

    query_names <- labelled_dat %>%
      filter(studytime == studytime_i, contrast == contrast_i,
             sigFDR == 1, !str_starts(metaboanalyst_label, "Total")) %>%
      distinct(metaboanalyst_label) %>%
      pull(metaboanalyst_label)

    if (!length(query_names)) next   # no FDR-significant compounds -> no cell

    cells[[length(cells) + 1]] <- list(
      study           = as.character(studytime_i),
      contrast        = contrast_i,
      direction       = "combined",   # pos + neg analysed together
      outcome_group   = target_group,
      query           = query_names,
      reference_names = reference_names
    )
  }
  cells
}

build_cells <- function(dat, outcome_group, measure = "ATE", pval_cutoff = 0.05) {
  # Note: `outcome_group`/`measure` name both the function arguments and the
  # data columns we filter on. Copy them to locals so the dplyr filter below
  # compares column-vs-value rather than column-vs-itself.
  target_group   <- outcome_group
  target_measure <- measure

  labelled_dat <- dat %>%
    filter(outcome_group == target_group, measure == target_measure) %>%
    # Exclude "Total *" roll-up features (composite quantities, not individual
    # metabolites) on the RAW label BEFORE the synonym map runs. Doing it here is
    # required because the map renames "total vitamin B1 (expressed as thiamin)"
    # -> "Thiamine", which would otherwise slip past the later
    # `!str_starts(metaboanalyst_label, "Total")` guard (and that guard is also
    # capitalization-fragile). Case-insensitive so "Total"/"total" behave alike.
    filter(!str_detect(label_f, regex("^\\s*total\\b", ignore_case = TRUE))) %>%
    mutate(metaboanalyst_label = apply_label_map(label_f))

  reference_names <- build_reference_names(labelled_dat)

  # One row per study-time x contrast combination present in this outcome group.
  study_contrast_grid <- labelled_dat %>%
    distinct(studytime, contrast) %>%
    arrange(studytime, contrast)

  cells <- list()
  for (i in seq_len(nrow(study_contrast_grid))) {
    studytime_i <- study_contrast_grid$studytime[i]
    contrast_i  <- study_contrast_grid$contrast[i]

    for (direction in names(DIRECTION_SIGNS)) {
      query_names <- build_query_names(
        labelled_dat, studytime_i, contrast_i,
        effect_sign = DIRECTION_SIGNS[[direction]], pval_cutoff = pval_cutoff)

      # Skip directions with no significant compounds, nothing to enrich.
      if (!length(query_names)) next

      cells[[length(cells) + 1]] <- list(
        study           = as.character(studytime_i),
        contrast        = contrast_i,
        direction       = direction,
        outcome_group   = target_group,
        query           = query_names,
        reference_names = reference_names
      )
    }
  }
  cells
}
