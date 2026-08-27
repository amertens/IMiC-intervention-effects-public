# =============================================================================
# src/metaboanalyst/tests/test-14-primary-reference-required.R
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
source("src/metaboanalyst/run-primary.R")

test_that("run_primary aborts loudly (does not silently degrade) when the reference file is missing", {
  old <- PRIMARY_CONFIG$reference_path
  on.exit(PRIMARY_CONFIG$reference_path <<- old, add = TRUE)
  PRIMARY_CONFIG$reference_path <<- "does/not/exist/reference.txt"
  err <- tryCatch(.primary_reference(), error = function(e) conditionMessage(e))
  expect_type(err, "character")
  expect_match(err, "reference metabolome not found", ignore.case = TRUE)
})

test_that("the real curated reference file resolves to a non-empty name vector", {
  ref <- .primary_reference()
  expect_true(is.character(ref))
  expect_gt(length(ref), 20)
})
