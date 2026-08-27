library(testthat)
# testthat::test_file() sets wd to the test folder; anchor to repo root.
repo_anchor <- function() {
  d <- getwd()
  for (i in 1:6) { if (dir.exists(file.path(d, "src/metaboanalyst"))) return(d); d <- dirname(d) }
  stop("could not locate repo root (src/metaboanalyst not found)")
}
setwd(repo_anchor())
source("src/metaboanalyst/build-supplementary-table.R")

# Uses the primary_combined outputs written by run-primary.R.
test_that("supplementary table consolidates cells with the right columns", {
  tab <- build_supplementary_table(write = FALSE)
  needed <- c("outcome_group", "arm_set", "study", "studytime", "direction", "contrast",
              "module", "pathway", "total", "hits", "expected", "raw_p",
              "fdr_native", "fdr_bh_pooled", "significant", "hit_features")
  expect_true(all(needed %in% names(tab)))
  expect_gt(nrow(tab), 0)
  expect_equal(unique(tab$outcome_group), "primary")
  expect_equal(unique(tab$arm_set), "combined")
})

test_that("Elicit-Up-1mo Nicotinate row is significant with 5 hit features", {
  tab <- build_supplementary_table(write = FALSE)
  row <- tab[tab$study == "Elicit" & tab$studytime == "1 mo." &
             tab$direction == "up" &
             tab$pathway == "Nicotinate and Nicotinamide Metabolism", ]
  expect_equal(nrow(row), 1)
  expect_equal(row$raw_p, 0.000367, tolerance = 1e-2)
  expect_true(row$significant)                       # fdr_native (0.036) < 0.05
  expect_equal(length(strsplit(row$hit_features, ";\\s*")[[1]]), 5)
})

test_that("significance uses native FDR, not pooled BH", {
  tab <- build_supplementary_table(write = FALSE)
  expect_equal(tab$significant, tab$fdr_native < 0.05)
})
