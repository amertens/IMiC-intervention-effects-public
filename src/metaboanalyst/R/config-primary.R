# config-primary.R, which cells to run for the primary outcome group.
# Adding a modality later = add a config like this one (module: "ora" or "pathway").
PRIMARY_CONFIG <- list(
  outcome_group = "primary",
  module        = "ora",     # Trenton's rule: ORA for primary
  measure       = "ATE",
  pval_cutoff   = 0.05,
  # Reference metabolome (background) for the ORA module. Trenton's curated file
  # holds the MetaboAnalyst-ID-conversion-MATCHED compound names (e.g.
  # "L-Tryptophan", "alpha-Tocopherol"), which match the SMPDB library far
  # better than the raw data labels do. Using it reproduces his reported
  # nicotinate structure (total 7 / hits 5); a data-derived name list matches
  # less completely (total 4 / hits 3). Auto-generating matched-name references
  # for other modalities is a later-phase task.
  # Rescued into the repo's own tree 2026-08-25 (was pointing into "trenton scripts/",
  # which has a paused, in-progress archival cleanup elsewhere in this repo -- see
  # Manuscript/CODE_AUDIT_2026-08-25.md #1). Byte-identical copy verified via diff.
  reference_path = "src/metaboanalyst/reference/Primary Outcomes Reference Metabolome.txt"
)
