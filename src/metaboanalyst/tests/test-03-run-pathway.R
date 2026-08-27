# =============================================================================
# src/metaboanalyst/tests/test-03-run-pathway.R
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
source("../R/run-pathway.R")

test_that("run_pathway reproduces the golden cell and exposes its workdir", {
  cmpd.vec <- c("Nicotinic acid mononucleotide", "Niacinamide", "Nudifloramide",
                "nicotinamide", "nicotinamide riboside",
                "nicotinamide adenine dinucleotide", "tryptophan")
  mSet <- run_pathway(cmpd.vec)                      # defaults: SMPDB, filter off
  expect_true(!is.null(mSet$imic_workdir))
  expect_true(file.exists(file.path(mSet$imic_workdir, "pathway_results.csv")))
  res <- read.csv(file.path(mSet$imic_workdir, "pathway_results.csv"),
                  row.names = 1, check.names = FALSE)
  row <- res["Nicotinate and Nicotinamide Metabolism", ]
  expect_equal(unname(row[["Raw p"]]), 1.2011e-05, tolerance = 1e-3)
  expect_equal(unname(row[["Hits"]]), 4)
})
