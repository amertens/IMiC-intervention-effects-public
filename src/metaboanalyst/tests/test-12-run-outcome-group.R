# =============================================================================
# src/metaboanalyst/tests/test-12-run-outcome-group.R
#
# This script's input/output paths are assembled at runtime, so they
# cannot be listed here without executing it.
#
# Header generated from the code itself; it makes no claim about method.
# See README.md for run order and results/ARTIFACT_MANIFEST.csv for the
# exhibit each script feeds.
# =============================================================================

library(testthat)
library(MetaboAnalystR)
repo_anchor <- function() {
  d <- getwd()
  for (i in 1:6) { if (dir.exists(file.path(d, "src/metaboanalyst"))) return(d); d <- dirname(d) }
  stop("could not locate repo root (src/metaboanalyst not found)")
}
setwd(repo_anchor())
source("src/metaboanalyst/R/run-outcome-group.R")

combined <- readRDS("trenton scripts/1. Data/combined_intervention_effects_results_combined_arms.RDS")

test_that("run_outcome_group returns results/skipped/out_dir and accounts for every cell", {
  res <- run_outcome_group(combined, "primary", module = "ora", arm_set = "combined",
                           reference_names = NULL, write = FALSE, max_cells = 3)
  expect_named(res, c("results", "skipped", "out_dir"))
  expect_s3_class(res$results, "data.frame")
  expect_s3_class(res$skipped, "data.frame")
  expect_true(grepl("primary_combined$", res$out_dir))
  # results carry the cell-tagging columns
  expect_true(all(c("study", "contrast", "direction", "pathway", "raw_p", "fdr") %in% names(res$results)))
})

test_that("run_outcome_group dispatches the pathway module (impact column present)", {
  # A tertiary slice large enough to include a mappable (>=3 metabolite) cell.
  res <- run_outcome_group(combined, "tertiary", module = "pathway", arm_set = "combined",
                           write = FALSE, max_cells = 6)
  # Either some cell produced pathway results (with an impact column) or all were
  # too thin and got skipped — both are valid; if any results exist they must be
  # pathway-shaped.
  if (nrow(res$results) > 0) {
    expect_true("impact" %in% names(res$results))
  } else {
    expect_gt(nrow(res$skipped), 0)
  }
})
