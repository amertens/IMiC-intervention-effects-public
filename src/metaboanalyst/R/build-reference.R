# build-reference.R, generate a MetaboAnalyst-matched reference metabolome from
# raw measured compound names. Automates Trenton's manual ID-conversion step
# (paste names into MetaboAnalyst's Compound ID Conversion, keep the "Match"
# column) so the ORA background is name-matched to the SMPDB/HMDB library.
#
# Returns the vector of matched compound names; optionally caches to a file.
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
