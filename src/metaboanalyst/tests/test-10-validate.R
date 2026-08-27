library(testthat)
library(MetaboAnalystR)
# testthat::test_file() sets wd to the test folder; anchor to repo root.
repo_anchor <- function() {
  d <- getwd()
  for (i in 1:6) { if (dir.exists(file.path(d, "src/metaboanalyst"))) return(d); d <- dirname(d) }
  stop("could not locate repo root (src/metaboanalyst not found)")
}
setwd(repo_anchor())
source("src/metaboanalyst/validate-against-trenton.R")

test_that("our golden pathway cell matches Trenton's downloaded CSV within 1%", {
  cmp <- validate_golden_pathway()
  expect_gt(nrow(cmp), 0)
  # every shared pathway's raw p within 1% relative tolerance of Trenton's
  expect_true(all(cmp$rel_diff < 0.01))
  # and the headline pathway is reproduced essentially exactly
  nic <- cmp[cmp$pathway == "Nicotinate and Nicotinamide Metabolism", ]
  expect_equal(nrow(nic), 1)
  expect_equal(nic$raw_p_ours, 1.2011e-05, tolerance = 1e-3)
})

test_that("our golden ORA cell matches Trenton's downloaded CSV", {
  cmp <- validate_golden_ora()
  expect_gt(nrow(cmp), 0)
  # Trenton's downloaded msea_ora_result.csv reports raw p to 3 significant
  # figures, so allow 2% relative tolerance (rounding-dominated) on every shared
  # pathway rather than the pathway gate's 1%.
  expect_true(all(cmp$rel_diff < 0.02))
  # The headline pathway (total 7 / hits 5, per config-primary.R) reproduces
  # Trenton's reported raw p of 0.000367.
  nic <- cmp[cmp$pathway == "Nicotinate and Nicotinamide Metabolism", ]
  expect_equal(nrow(nic), 1)
  expect_equal(nic$raw_p_ours, 0.000367, tolerance = 2e-2)
})
