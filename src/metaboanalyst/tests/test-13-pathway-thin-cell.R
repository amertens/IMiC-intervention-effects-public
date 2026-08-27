# =============================================================================
# src/metaboanalyst/tests/test-13-pathway-thin-cell.R
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

test_that("run_pathway raises a clean catchable error on a too-few-metabolites cell", {
  # Even with current.msg seeded (as an earlier run_ora would leave it), the
  # too-few path must NOT segfault -- it must raise a catchable, informative error.
  assign("current.msg", character(0), envir = .GlobalEnv)   # simulate prior run_ora seeding
  err <- tryCatch(run_pathway(c("Cholesterol", "Glucose")), error = function(e) conditionMessage(e))
  expect_type(err, "character")
  expect_match(err, "too few mappable metabolites", ignore.case = TRUE)
})

test_that("golden pathway cell still reproduces after the robustness fix", {
  g <- c("Nicotinic acid mononucleotide", "Niacinamide", "Nudifloramide", "nicotinamide",
         "nicotinamide riboside", "nicotinamide adenine dinucleotide", "tryptophan")
  res <- harvest_results(run_pathway(g, "smpdb"))
  nic <- res[grepl("Nicotinate", res$pathway), ][1, ]
  expect_equal(nic$raw_p, 1.2011e-05, tolerance = 1e-3)
})
