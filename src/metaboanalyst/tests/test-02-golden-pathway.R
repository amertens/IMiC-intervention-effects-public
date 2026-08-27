# Reproduces the golden Elicit-Up-1mo pathway cell verbatim from its downloaded Rhistory.R.
# Library confirmed as SMPDB (SetSMPDB.PathLib) from the downloaded Rhistory + SMPDB-style pathway names.
library(testthat)
library(MetaboAnalystR)

test_that("golden Elicit-Up-1mo pathway cell matches downloaded pathway_results.csv", {
  cmpd.vec <- c("Nicotinic acid mononucleotide", "Niacinamide", "Nudifloramide",
                "nicotinamide", "nicotinamide riboside",
                "nicotinamide adenine dinucleotide", "tryptophan")

  workdir <- tempfile("golden_"); dir.create(workdir); old <- setwd(workdir); on.exit(setwd(old))

  # MetaboAnalystR 4.3.0's InitDataObjects() has a self-referential default
  # argument (default.dpi = default.dpi), which throws "promise already under
  # evaluation" if the 4th arg is omitted. Passing an explicit dpi (unrelated
  # to pathway-analysis numbers) works around this without touching the
  # compound list, library, or filter setting.
  mSet <- InitDataObjects("conc", "pathora", FALSE, 72)
  mSet <- Setup.MapData(mSet, cmpd.vec)
  mSet <- CrossReferencing(mSet, "name")
  mSet <- CreateMappingResultTable(mSet)
  mSet <- SetSMPDB.PathLib(mSet, "hsa")
  mSet <- SetOrganism(mSet, "hsa")
  mSet <- SetMetabolomeFilter(mSet, FALSE)   # matches the golden Rhistory (filter OFF)
  mSet <- CalculateOraScore(mSet, "rbc", "hyperg")

  # mSet$analSet$ora.mat keys rows by internal SMPDB IDs (e.g. "SMP00048"), not
  # pathway names. CalculateOraScore() also writes "pathway_results.csv" to the
  # working directory with rownames remapped to full pathway names via
  # current.kegglib$path.ids — that file is the same artifact Trenton
  # downloaded, so read it back instead of the raw ID-keyed matrix.
  res <- read.csv(file.path(workdir, "pathway_results.csv"), row.names = 1,
                   check.names = FALSE)
  expect_true("Nicotinate and Nicotinamide Metabolism" %in% rownames(res))
  row <- res["Nicotinate and Nicotinamide Metabolism", ]
  # Column name for raw p may be "Raw p" or "Raw.p" after as.data.frame — handle both.
  pcol <- intersect(c("Raw p","Raw.p"), colnames(res))[1]
  fcol <- intersect(c("FDR","FDR.q"), colnames(res))[1]
  expect_equal(unname(row[[pcol]]), 1.2011e-05, tolerance = 1e-3)
  expect_equal(unname(row[[fcol]]), 1.1891e-03, tolerance = 1e-3)
  expect_equal(unname(row[["Hits"]]),  4)
  expect_equal(unname(row[["Total"]]), 32)
})
