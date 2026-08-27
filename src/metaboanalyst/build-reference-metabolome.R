# build-reference-metabolome.R -- reproducibly derive the reference metabolome
# used as the enrichment background for the manuscript's tertiary (Fig 5B) and
# optional primary (Fig 3B) MSEA panels.
#
# NOTE ON NAMING: this file was previously (mis)named "whole-metabolome". It is
# NOT the whole human metabolome -- it is the author's saved reference metabolome
# (metabolites + lipids). The larger submission-matching variant is the QER
# reference metabolome (src/metaboanalyst/reference/refMetabolomeForQER.csv, 1,268
# names), which reproduces the submitted Fig 5B; this 1,593-HMDB variant is the
# smaller relative used by run-tertiary-msea.R's default background.
#
# Source of truth: the author's reference_metabolome.xlsx, a single column of
# 1,593 HMDB accession IDs (1,273 unique). MetaboAnalyst's metabolite-set library
# is keyed by compound NAME, so we cross-reference the HMDB IDs to library names;
# 615 of the 1,273 resolve to a named compound the SMPDB library recognises. The
# unmatched 658 are HMDB entries with no member in any named metabolite set -- they
# cannot contribute to name-based enrichment, so the effective background is these
# 615 names.
#
# Inputs (committed):  src/metaboanalyst/reference/reference_metabolome_1593_hmdb.txt
# Output (committed):   src/metaboanalyst/reference/reference_metabolome_1593_matched_names.txt
#
# Run from repo root:  Rscript src/metaboanalyst/build-reference-metabolome.R
suppressMessages({ library(MetaboAnalystR) })

HMDB_IN   <- "src/metaboanalyst/reference/reference_metabolome_1593_hmdb.txt"
NAMES_OUT <- "src/metaboanalyst/reference/reference_metabolome_1593_matched_names.txt"

build_reference_metabolome <- function(hmdb_in = HMDB_IN, names_out = NAMES_OUT) {
  ids <- unique(trimws(readLines(hmdb_in, warn = FALSE)))
  ids <- ids[ids != ""]
  message("input HMDB IDs: ", length(ids))

  # AddErrMsg (unmatched-name reporting) reads these globals off the public web.
  for (g in c("current.msg", "err.vec"))
    if (!exists(g, envir = .GlobalEnv)) assign(g, character(0), envir = .GlobalEnv)

  # CrossReferencing writes name_map.csv to the wd; isolate + restore.
  wd <- tempfile("hmdbmap_"); dir.create(wd); ow <- setwd(wd); on.exit(setwd(ow), add = TRUE)
  mSet <- InitDataObjects("conc", "msetora", FALSE, 150)
  mSet <- Setup.MapData(mSet, ids)
  mSet <- CrossReferencing(mSet, "hmdb")          # input type = HMDB accession IDs
  mSet <- CreateMappingResultTable(mSet)
  nm <- utils::read.csv(file.path(wd, "name_map.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  setwd(ow)

  matched <- unique(nm$Match[!is.na(nm$Match) & nm$Match != ""])
  message("matched library NAMES: ", length(matched), " of ", length(ids))

  dir.create(dirname(names_out), recursive = TRUE, showWarnings = FALSE)
  writeLines(matched, names_out)
  message("wrote ", names_out)
  invisible(matched)
}

if (sys.nframe() == 0) invisible(build_reference_metabolome())
