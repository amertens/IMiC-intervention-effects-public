# =============================================================================
# src/metaboanalyst/tests/test-07-build-cells.R
#
# This script's input/output paths are assembled at runtime, so they
# cannot be listed here without executing it.
#
# Header generated from the code itself; it makes no claim about method.
# See README.md for run order and results/ARTIFACT_MANIFEST.csv for the
# exhibit each script feeds.
# =============================================================================

library(testthat)
repo_anchor <- function() {
  d <- getwd()
  for (i in 1:6) { if (dir.exists(file.path(d, "src/metaboanalyst"))) return(d); d <- dirname(d) }
  stop("could not locate repo root (src/metaboanalyst not found)")
}
setwd(repo_anchor())
source("src/metaboanalyst/R/build-cells.R")

combined <- readRDS("trenton scripts/1. Data/combined_intervention_effects_results_combined_arms.RDS")

test_that("build_cells recovers the known Elicit-Up-1mo query list", {
  cells <- build_cells(combined, outcome_group = "primary")
  idx <- which(vapply(cells, function(c)
    c$study == "Elicit (1 mo.)" && c$direction == "up" && c$contrast == "Nico", logical(1)))
  expect_length(idx, 1)
  cell <- cells[[idx]]
  expect_setequal(
    cell$query,
    c("Niacinamide", "nicotinamide", "nicotinamide adenine dinucleotide",
      "nicotinamide riboside", "Nicotinic acid mononucleotide", "Nudifloramide", "tryptophan")
  )
})

test_that("reference metabolome excludes Total* and is non-empty", {
  cells <- build_cells(combined, outcome_group = "primary")
  expect_false(any(grepl("^Total", cells[[1]]$reference_names)))
  expect_gt(length(cells[[1]]$reference_names), 20)
})
