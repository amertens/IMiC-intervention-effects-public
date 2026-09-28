# =============================================================================
# harvest.R
#
# Helpers that turn a MetaboAnalystR mSet from run_ora() or run_pathway() into
# tidy tables: harvest_results() (one row per pathway: total, hits, expected,
# raw P, FDR, and impact for pathway runs), harvest_membership() (pathway ->
# hit compounds), and apply_pathway_size_floor() (reporting filter on pathway
# size). Sourced by run-tertiary-msea-dual.R and run-untargeted-msea.R.
#
# Inputs : the results CSV MetaboAnalystR writes into the run's working directory
# Outputs: none (returns tibbles)
# =============================================================================
suppressMessages({ library(dplyr); library(tibble) })

# apply_pathway_size_floor: drop pathways with fewer than `min_size` members
# (in the reference-restricted library) from the reported table. 1-2 member sets
# give enrichment ratios in the hundreds or thousands that are noise.
#
# MetaboAnalystR's own per-cell BH FDR (`fdr_native`) is kept unchanged. It is
# computed over the full reference-restricted SMPDB library (about 89-99 sets per
# cell), which is the right multiple-testing family. BH is not recomputed over
# the surviving rows: the harvested table holds only sets with >= 1 hit (about 26
# of 89), so re-running p.adjust() on it would shrink the family and inflate
# significance. The floor is therefore a reporting filter, and `fdr_native` is
# slightly conservative because the dropped small sets still count toward m.
# Used by the tertiary (Fig 5B; min_size = 1) and untargeted (Fig 6A; min_size = 3) runs.
#
# `group_cols` is unused; it is kept so existing calls keep working.
apply_pathway_size_floor <- function(tab, group_cols = NULL, min_size = 3) {
  dplyr::filter(tab, .data$total >= min_size)
}

# Locate and read the single results CSV that a run_ora/run_pathway mSet wrote to
# its workdir. Both harvest_results and harvest_membership read the SAME file, so
# the search pattern lives here once to keep them aligned.
read_imic_results_csv <- function(mSet) {
  stopifnot(!is.null(mSet$imic_workdir))
  csv <- list.files(mSet$imic_workdir,
                    pattern = "pathway_results\\.csv|msea_ora_result\\.csv|ora\\.csv",
                    full.names = TRUE)
  if (!length(csv)) stop("no results CSV found in imic_workdir")
  read.csv(csv[1], row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
}

harvest_results <- function(mSet) {
  results_table <- read_imic_results_csv(mSet)
  present_cols  <- colnames(results_table)

  # The two modules label the same quantities with different column headers, so
  # we look up each field by trying a list of candidate names and taking the
  # first that exists. Returns NA when none is present.
  #  - pathora's pathway_results.csv uses "Total"/"Hits"/"Expected" (title case);
  #  - msetora's msea_ora_result.csv uses "total"/"hits"/"expected" (lowercase).
  first_present_col <- function(candidates) {
    found <- candidates[candidates %in% present_cols]
    if (length(found)) found[1] else NA_character_
  }
  total_col    <- first_present_col(c("Total", "total"))
  hits_col     <- first_present_col(c("Hits", "hits"))
  expected_col <- first_present_col(c("Expected", "expected"))
  raw_p_col    <- first_present_col(c("Raw p", "Raw.p", "RawP"))
  # Accept only BH/FDR column names. Do not fall back to a Holm column: Holm is a
  # different (FWER) procedure, and storing a Holm-adjusted value as `fdr_native`
  # would silently change the significance definition downstream (significant =
  # fdr_native < 0.05). If no FDR column is present, `fdr` stays NA (surfaced), not
  # a mislabeled Holm value.
  fdr_col      <- first_present_col(c("FDR", "FDR.q"))

  # Pull a numeric column by header, tolerating absence (NA) and non-numeric
  # text (suppressed coercion warnings -> NA).
  numeric_col <- function(col_name) {
    if (is.na(col_name)) return(NA_real_)
    suppressWarnings(as.numeric(results_table[[col_name]]))
  }

  out <- tibble(
    pathway  = rownames(results_table),
    total    = numeric_col(total_col),
    hits     = numeric_col(hits_col),
    expected = numeric_col(expected_col),
    raw_p    = numeric_col(raw_p_col),
    fdr      = numeric_col(fdr_col)
  )
  # Impact is a pathway-topology metric present only in the pathora output.
  if ("Impact" %in% present_cols) {
    out$impact <- suppressWarnings(as.numeric(results_table[["Impact"]]))
  }
  arrange(out, raw_p)
}

# harvest_membership: pathway -> hit-feature membership table.
#
# mSet$analSet$ora.hits is a list, one entry per pathway in the full library
# (including pathways later dropped from the results table, e.g. by min/max-size
# filtering). Pathways with no hits have a length-0 vector.
#
# The keying and vector shape differ by module:
#  - pathora (run_pathway):  list keyed by internal SMPDB IDs (e.g. "SMP00048");
#    each entry is a NAMED character vector: names(...) = hit feature names,
#    values = their HMDB IDs.
#  - msetora (run_ora):      list keyed directly by pathway NAME (matches
#    rownames(ora.mat)/the results CSV already); each entry is an UNNAMED
#    character vector whose VALUES are the hit feature names.
# We handle both: prefer names(x) when present (pathora), else fall back to the
# vector's own values (msetora).
#
# The name-keyed results CSV (pathway_results.csv) and the ID-keyed ora.mat share
# row order (both built from the same filtered/ordered pathway set), so we zip
# them together to get an ID -> full pathway name lookup, then use that to keep
# only the hits belonging to pathways that made it into the final results table.
harvest_membership <- function(mSet) {
  results_table <- read_imic_results_csv(mSet)

  ora.mat <- mSet$analSet$ora.mat
  hits    <- mSet$analSet$ora.hits
  if (is.null(ora.mat)) stop("mSet$analSet$ora.mat not found")
  if (is.null(hits))    stop("mSet$analSet$ora.hits not found")
  if (nrow(results_table) != nrow(ora.mat)) {
    stop("pathway_results.csv and ora.mat row counts differ; cannot align ID->name mapping")
  }
  # results_table is keyed by full pathway name, ora.mat by internal ID, and the
  # two share row order, so zipping their rownames gives an ID -> name lookup.
  pathway_id_to_name <- setNames(rownames(results_table), rownames(ora.mat))

  membership_rows <- lapply(names(hits), function(pathway_id) {
    # Single-bracket lookup yields NA (not an error) for IDs that were filtered
    # out of the final results table; drop those.
    pathway_name <- unname(pathway_id_to_name[pathway_id])
    if (is.na(pathway_name)) return(NULL)

    # pathora entries are named (names = hit features); msetora entries are
    # unnamed, so the values themselves are the hit feature names.
    feature_names <- names(hits[[pathway_id]])
    if (is.null(feature_names)) feature_names <- unname(hits[[pathway_id]])
    if (length(feature_names) == 0) return(NULL)  # pathway had no hits

    tibble(pathway = pathway_name, feature = feature_names)
  })
  bind_rows(membership_rows)
}
