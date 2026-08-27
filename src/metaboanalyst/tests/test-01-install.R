# =============================================================================
# src/metaboanalyst/tests/test-01-install.R
#
# This script's input/output paths are assembled at runtime, so they
# cannot be listed here without executing it.
#
# Header generated from the code itself; it makes no claim about method.
# See README.md for run order and results/ARTIFACT_MANIFEST.csv for the
# exhibit each script feeds.
# =============================================================================

library(testthat)
test_that("MetaboAnalystR loads and core functions exist", {
  expect_true(requireNamespace("MetaboAnalystR", quietly = TRUE))
  library(MetaboAnalystR)
  for (fn in c("InitDataObjects", "Setup.MapData", "CrossReferencing",
               "CreateMappingResultTable", "CalculateOraScore",
               "SetKEGG.PathLib", "SetMetabolomeFilter", "CalculateHyperScore")) {
    expect_true(exists(fn), info = fn)
  }
})
