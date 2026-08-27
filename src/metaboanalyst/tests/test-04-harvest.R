# =============================================================================
# src/metaboanalyst/tests/test-04-harvest.R
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
source("src/metaboanalyst/R/run-pathway.R")
source("src/metaboanalyst/R/harvest.R")

test_that("harvest_results returns a tidy results tibble", {
  cmpd.vec <- c("Nicotinic acid mononucleotide", "Niacinamide", "Nudifloramide",
                "nicotinamide", "nicotinamide riboside",
                "nicotinamide adenine dinucleotide", "tryptophan")
  mSet <- run_pathway(cmpd.vec)
  results <- harvest_results(mSet)
  expect_true(all(c("pathway", "total", "hits", "raw_p", "fdr") %in% names(results)))
  top <- results[results$pathway == "Nicotinate and Nicotinamide Metabolism", ]
  expect_equal(top$raw_p, 1.2011e-05, tolerance = 1e-3)
  # results should be sorted ascending by raw_p (most significant first)
  expect_equal(results$pathway[1], "Nicotinate and Nicotinamide Metabolism")
})

test_that("harvest_membership maps pathways to their hit features", {
  cmpd.vec <- c("Nicotinic acid mononucleotide", "Niacinamide", "Nudifloramide",
                "nicotinamide", "nicotinamide riboside",
                "nicotinamide adenine dinucleotide", "tryptophan")
  mSet <- run_pathway(cmpd.vec)
  membership <- harvest_membership(mSet)
  expect_true(all(c("pathway", "feature") %in% names(membership)))
  nic <- membership[membership$pathway == "Nicotinate and Nicotinamide Metabolism", ]
  expect_gt(nrow(nic), 0)   # the nicotinate pathway must list its hit features
})
