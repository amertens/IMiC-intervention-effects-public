# build-supplementary-table.R, consolidate one outcome-group x arm-set's cells
# into a tidy supplementary table.
#
# Works for any group directory written by run_outcome_group()
# (results/metaboanalyst/<group>_<arm>/<cell>/{results.csv,membership.csv}).
# Significance uses MetaboAnalyst's own per-analysis FDR (fdr_native), the
# package/BH correction Trenton used. fdr_bh_pooled (BH across all cells pooled)
# is informational only and does NOT drive the flag. Pathway-module tables carry
# an extra `impact` column; ORA tables do not.
suppressMessages({ library(dplyr); library(readr); library(stringr); library(tidyr); library(purrr) })

DEFAULT_GROUP_DIR <- "results/metaboanalyst/primary_combined"

# Pathway-size floor. Set to 1 (NO floor) to match Trenton's ORA figures and Fig 5B's
# dual runner, which apply no minimum set-size filter -- so the reported ORA tables and
# Fig S8 retain the same 1-2 member sets his submission shows. (Empirically inert on the
# current data: 0 FDR-significant and 0 nominally-significant primary/tertiary pathways
# have total<3, so this changes no significant point/label; it only removes a silent
# divergence from the submission and the earlier 5B-vs-S8/6A inconsistency.) 1-2 member
# sets can carry noisy enrichment ratios (ER = hits/expected in the hundreds), so they
# should be interpreted cautiously where present. `fdr_native` is left unchanged.
MIN_PATHWAY_SIZE <- 1

# Cell labels arrive as "Study (collection-time window)", e.g. "Elicit (1 month)".
# These two helpers split that label: .study_name strips the parenthetical to the
# bare study name; .study_time returns only the text inside the parentheses.
.study_name <- function(label) str_trim(str_replace(label, "\\s*\\(.*$", ""))
.study_time <- function(label) {
  inside_parens <- str_replace(label, "^[^(]*\\(", "")   # drop everything up to "("
  str_trim(str_replace(inside_parens, "\\)\\s*$", ""))   # drop the closing ")"
}

build_supplementary_table <- function(group_dir = DEFAULT_GROUP_DIR, write = TRUE) {
  # The directory is named "<outcome_group>_<arm_set>"; arm_set is the final
  # underscore-delimited token and outcome_group is everything before it.
  base <- basename(group_dir)
  arm_set       <- str_extract(base, "[^_]+$")
  outcome_group <- str_replace(base, "_[^_]+$", "")

  cell_dirs <- list.dirs(group_dir, recursive = FALSE)
  cell_dirs <- cell_dirs[file.exists(file.path(cell_dirs, "results.csv"))]
  if (!length(cell_dirs)) {
    message("build_supplementary_table: no cell results in ", group_dir,
            " (all cells skipped) -- no table written.")
    return(invisible(NULL))
  }

  results <- map_dfr(cell_dirs, ~ read_csv(file.path(.x, "results.csv"), show_col_types = FALSE))
  membership <- map_dfr(cell_dirs, ~ read_csv(file.path(.x, "membership.csv"), show_col_types = FALSE))
  # Pathway-module runs carry an `impact` column; ORA runs do not. Presence of
  # that column tells us which module produced these cells.
  module <- if ("impact" %in% names(results)) "pathway" else "ora"

  # For each cell x pathway, collapse the enriched ("hit") compounds into one
  # semicolon-joined string so the table can show what drove each enrichment.
  hit_features_per_pathway <- membership %>%
    group_by(study, contrast, direction, pathway) %>%
    summarise(hit_features = paste(feature, collapse = "; "), .groups = "drop")

  supp_table <- results %>%
    rename(fdr_native = fdr) %>%
    # Size floor (see MIN_PATHWAY_SIZE): drop 1-2 member noise sets from the
    # reported table and the `significant` flag, matching the MSEA convention.
    filter(total >= MIN_PATHWAY_SIZE) %>%
    left_join(hit_features_per_pathway, by = c("study", "contrast", "direction", "pathway")) %>%
    mutate(
      hit_features  = replace_na(hit_features, ""),
      outcome_group = outcome_group,
      arm_set       = arm_set,
      module        = module,
      # `study` still holds the full "Name (time window)" label at this point, so
      # derive studytime BEFORE overwriting `study` with the bare name below.
      studytime     = .study_time(study),
      fdr_bh_pooled = p.adjust(raw_p, method = "BH"),
      significant   = fdr_native < 0.05,
      study         = .study_name(study)
    )

  base_cols <- c("outcome_group", "arm_set", "study", "studytime", "direction",
                 "contrast", "module", "pathway", "total", "hits", "expected",
                 "raw_p", "fdr_native", "fdr_bh_pooled", "significant")
  impact_col <- if ("impact" %in% names(supp_table)) "impact" else character(0)
  supp_table <- supp_table %>%
    select(all_of(base_cols), all_of(impact_col), hit_features) %>%
    arrange(study, studytime, direction, raw_p)

  if (isTRUE(write)) {
    prefix <- file.path(group_dir, paste0(outcome_group, "_", arm_set))
    write_csv(supp_table, paste0(prefix, "_supplementary_table.csv"))
    write_csv(filter(supp_table, significant), paste0(prefix, "_significant.csv"))
    writeLines(.table_readme(outcome_group, arm_set, module), file.path(group_dir, "TABLE-README.md"))
  }
  supp_table
}

.table_readme <- function(outcome_group, arm_set, module) {
  c(
    sprintf("# %s / %s enrichment supplementary table (%s module)", outcome_group, arm_set, module),
    "",
    "`<group>_<arm>_supplementary_table.csv` holds one row per (analysis cell x pathway).",
    "`<group>_<arm>_significant.csv` is the `significant == TRUE` subset.",
    "",
    "## Columns",
    "- outcome_group, arm_set, study, studytime, direction, contrast, module: identify the cell.",
    "- pathway: SMPDB metabolite-set / pathway name.",
    "- total: pathway members in the background. hits: query compounds hitting it.",
    "- expected: hits expected by chance. raw_p: over-representation p-value.",
    "- **fdr_native**: MetaboAnalyst's own per-analysis FDR (Benjamini-Hochberg within the",
    "  analysis): the correction used for significance, matching Trenton's approach.",
    "- fdr_bh_pooled: BH across ALL pathways x cells pooled, informational only.",
    "- **significant**: fdr_native < 0.05.",
    sprintf("- Pathways with < %d members in the (reference-restricted) library are dropped", MIN_PATHWAY_SIZE),
    "  as noise (reporting filter applied after p-values; fdr_native left unchanged).",
    if (module == "pathway") "- impact: pathway topology impact (relative-betweenness).",
    "- hit_features: enriched compounds for that (cell, pathway), semicolon-joined.",
    "",
    "Regenerate: `Rscript src/metaboanalyst/build-supplementary-table.R <group_dir>`"
  )
}

if (sys.nframe() == 0) {
  a <- commandArgs(trailingOnly = TRUE)
  invisible(build_supplementary_table(if (length(a)) a[1] else DEFAULT_GROUP_DIR, write = TRUE))
}
