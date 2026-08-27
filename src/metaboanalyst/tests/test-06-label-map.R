library(testthat)
# testthat::test_file() sets wd to the test file's folder; anchor to repo root
# so repo-root-relative source()/readRDS() paths resolve however this is launched.
repo_anchor <- function() {
  d <- getwd()
  for (i in 1:6) { if (dir.exists(file.path(d, "src/metaboanalyst"))) return(d); d <- dirname(d) }
  stop("could not locate repo root (src/metaboanalyst not found)")
}
setwd(repo_anchor())
source("src/metaboanalyst/R/label-map.R")

test_that("apply_label_map renames known synonyms and leaves others untouched", {
  x <- c("vitamin B3 ( expressed as NAM)", "vitamin B2 (expressed as riboflavin)",
         "a.tocopherol", "free thiamin", "potassium")
  out <- apply_label_map(x)
  expect_equal(out, c("Niacinamide", "Riboflavin", "Alpha-tocopherol", "Thiamine", "potassium"))
})

test_that("apply_label_map handles the tocopherol and vitamin A anchored patterns", {
  expect_equal(apply_label_map("g.tocopherol"), "Gamma-tocopherol")
  expect_equal(apply_label_map("vitamin.a"), "Vitamin A")
  expect_equal(apply_label_map("nicotinamide mononucleotide"), "Nicotinic acid mononucleotide")
})
