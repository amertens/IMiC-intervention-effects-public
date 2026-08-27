# =============================================================================
# src/metaboanalyst/tests/test-08-run-ora.R
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
source("src/metaboanalyst/R/run-ora.R")
source("src/metaboanalyst/R/harvest.R")

GOLDEN7 <- c("Nicotinic acid mononucleotide","Niacinamide","Nudifloramide","nicotinamide",
             "nicotinamide riboside","nicotinamide adenine dinucleotide","tryptophan")

test_that("run_ora (no reference) runs and flags the nicotinate set", {
  mSet <- run_ora(GOLDEN7, reference_names = NULL)
  res  <- harvest_results(mSet)
  expect_true(all(c("pathway","raw_p","fdr") %in% names(res)))
  expect_true(any(grepl("Nicotinate", res$pathway)))
  nic <- res[grepl("Nicotinate", res$pathway), ][1, ]
  expect_equal(nic$hits, 5)
})

test_that("run_ora with a reference metabolome restricts the background", {
  mSet <- run_ora(GOLDEN7, reference_names = GOLDEN7)   # tiny reference just to exercise the path
  res  <- harvest_results(mSet)
  expect_true(any(grepl("Nicotinate", res$pathway)))
})

test_that("harvest_membership works on an ORA mSet", {
  mSet <- run_ora(GOLDEN7, reference_names = NULL)
  mem  <- harvest_membership(mSet)
  expect_true(all(c("pathway","feature") %in% names(mem)))
  expect_gt(nrow(mem[grepl("Nicotinate", mem$pathway), ]), 0)
})
