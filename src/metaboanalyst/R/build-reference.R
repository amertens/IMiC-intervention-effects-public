# =============================================================================
# build-reference.R
#
# Helper: build_matched_reference() turns raw measured compound names into
# MetaboAnalyst-matched names, the scripted form of the web tool's Compound ID
# Conversion step (keep the "Match" column), so an ORA background is
# name-matched to the SMPDB/HMDB library. Returns the matched names and can
# cache them to a file. Sourced by run-tertiary-msea-dual.R; the published runs
# use the saved reference lists in src/metaboanalyst/reference/ instead.
#
# Inputs : none (takes a character vector)
# Outputs: optional cache file given by `cache_path`
# =============================================================================
suppressMessages({ library(MetaboAnalystR) })

build_matched_reference <- function(raw_names, cache_path = NULL, refresh = FALSE, lipid = FALSE) {
  if (!refresh && !is.null(cache_path) && file.exists(cache_path)) {
    return(trimws(readLines(cache_path, warn = FALSE)))
  }
  # AddErrMsg (called to report unmatched names) reads these globals off-web.
  for (global_name in c("current.msg", "err.vec")) {
    if (!exists(global_name, envir = .GlobalEnv)) assign(global_name, character(0), envir = .GlobalEnv)
  }
  # CrossReferencing writes name_map.csv into the working directory, so run in a
  # throwaway dir and restore the caller's wd on exit.
  workdir <- tempfile("mapref_"); dir.create(workdir)
  original_wd <- setwd(workdir); on.exit(setwd(original_wd), add = TRUE)

  mSet <- InitDataObjects("conc", "msetora", FALSE, 150)
  mSet <- Setup.MapData(mSet, unique(raw_names))
  mSet <- CrossReferencing(mSet, "name", lipid = lipid)
  mSet <- CreateMappingResultTable(mSet)

  # name_map.csv has one row per input name; its "Match" column holds the
  # library-recognized name (blank/NA when the name did not map). Keep the
  # distinct non-blank matches as the reference metabolome.
  name_map <- utils::read.csv(file.path(workdir, "name_map.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  matched <- unique(name_map$Match[!is.na(name_map$Match) & name_map$Match != ""])

  if (!is.null(cache_path)) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(matched, cache_path)
  }
  matched
}
