# =============================================================================
# build-reference-metabolome.R
#
# Documents how the 615-name reference metabolome was built. That list is the
# enrichment background for the untargeted milk ORA (run-untargeted-msea.R;
# Fig 6A, Table S5). The source is the study reference metabolome: 1,593 HMDB
# accession IDs (1,273 unique). MetaboAnalyst's metabolite-set library is keyed
# by compound name, so the HMDB IDs are cross-referenced to library names; 615
# resolve to a named compound. The other 658 have no member in any named
# metabolite set and cannot contribute to name-based enrichment. (The tertiary
# ORA uses a different, 1,268-name reference: reference/refMetabolomeForQER.csv.)
#
# Inputs : src/metaboanalyst/reference/reference_metabolome_1593_hmdb.txt
# Outputs: src/metaboanalyst/reference/reference_metabolome_1593_matched_names.txt
#
# MetaboAnalystR downloads its compound library from metaboanalyst.ca at run time.
# Run from repo root: Rscript src/metaboanalyst/build-reference-metabolome.R
# =============================================================================
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
