# =============================================================================
# src/metaboanalyst/tests/test-09-run-primary.R
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

test_that("run_primary produces a tidy results tibble for the primary cells", {
  out <- run_primary(write = FALSE, max_cells = 3)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("study","contrast","direction","pathway","raw_p","fdr") %in% names(out)))
  expect_gt(nrow(out), 0)
})
