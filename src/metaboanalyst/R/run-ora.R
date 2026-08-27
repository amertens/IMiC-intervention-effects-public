# run-ora.R — thin wrapper around MetaboAnalystR's msetora (over-representation)
# engine. Mirrors run-pathway.R's structure. Reproduces the sequence Task 5
# verified against Trenton's golden 7-compound cell, including the required
# session patch for CalculateHyperScore()'s off-public-web reporting-side-effect
# crash (see patch_msetora_export_bug() below).

# ---------------------------------------------------------------------------
# Known off-public-web bug in CalculateHyperScore() (msetora path only): it
# ends with ExportOraMembershipJson() and PlotORAMembership(NA, ...), both of
# which omit mSetObj. .get.mSet(NA) short-circuits to the NA literal off the
# public web (instead of falling back to a global mSetObj), so the very next
# line in each helper (`mSetObj$analSet$ora.mat`) throws "$ operator is
# invalid for atomic vectors". This happens unconditionally, even in the
# plain filter-OFF/no-reference case. Both helpers are reporting/plotting
# side effects only (a JSON membership dump and a PNG heatmap) that run AFTER
# ora.mat is already computed and stored -- no-op'ing them is numerically
# harmless (empirically verified in Task 5).
#
# AddErrMsg() also reads an uninitialized `current.msg` global off the public
# web. We do NOT seed it: seeding lets the too-few-metabolites error path
# proceed into C code that SEGFAULTS the process (same failure fixed in
# run_pathway). Instead we un-seed current.msg before CalculateHyperScore so the
# error path raises a catchable R error, which we convert to a clean skip reason.
# ---------------------------------------------------------------------------

# Replace the two crashing reporting helpers with harmless no-ops in the
# MetaboAnalystR namespace (see the note above for why this is safe). We must
# unlock each binding before reassigning it and re-lock it afterwards, because
# namespace bindings are locked by default.
patch_msetora_export_bug <- function() {
  metabo_ns <- asNamespace("MetaboAnalystR")
  for (fn_name in c("ExportOraMembershipJson", "PlotORAMembership")) {
    unlockBinding(fn_name, metabo_ns)
    assign(fn_name, function(mSetObj = NA, ...) invisible(NULL), envir = metabo_ns)
    lockBinding(fn_name, metabo_ns)
  }
  invisible(TRUE)
}

run_ora <- function(cmpd.vec,
                     reference_names = NULL,
                     mset_lib = "smpdb_pathway",
                     lipid = FALSE,   # TRUE = web "lipids" feature type (lipid_compound_db name match)
                     dpi = 150) {
  patch_msetora_export_bug()
  stopifnot(length(cmpd.vec) >= 1)

  # Seed current.msg / err.vec for the SETUP stages: Setup.HMDBReferenceMetabolome
  # (and other pre-score steps) call AddErrMsg() to report unmatched reference
  # names, which reads both globals. We remove current.msg again just before
  # scoring (segfault guard); err.vec is harmless to leave seeded.
  for (global_name in c("current.msg", "err.vec")) {
    if (!exists(global_name, envir = .GlobalEnv)) assign(global_name, character(0), envir = .GlobalEnv)
  }

  # MetaboAnalystR writes result CSVs to the working directory, so give each run
  # its own throwaway dir and restore the caller's wd on exit.
  workdir <- tempfile("ora_"); dir.create(workdir)
  original_wd <- setwd(workdir); on.exit(setwd(original_wd), add = TRUE)

  # dpi arg is mandatory in MetaboAnalystR 4.3.0 (self-referential default bug).
  mSet <- InitDataObjects("conc", "msetora", FALSE, dpi)
  mSet <- Setup.MapData(mSet, cmpd.vec)
  # lipid = TRUE swaps compound_db.qs -> lipid_compound_db.qs (the web "lipids"
  # feature type). Feed LIPID MAPS-formatted names (see R/lipid-name-map.R).
  mSet <- CrossReferencing(mSet, "name", lipid = lipid)
  mSet <- CreateMappingResultTable(mSet)

  # With a reference metabolome, restrict the enrichment background to those
  # compounds (filter ON); without one, use the full library (filter OFF).
  if (!is.null(reference_names)) {
    reference_path <- file.path(workdir, "reference_metabolome.txt")
    writeLines(reference_names, reference_path)
    mSet <- Setup.HMDBReferenceMetabolome(mSet, reference_path)
    mSet <- SetMetabolomeFilter(mSet, TRUE)
  } else {
    mSet <- SetMetabolomeFilter(mSet, FALSE)
  }

  mSet <- SetCurrentMsetLib(mSet, mset_lib, 2)

  # Robust too-few-metabolites handling (mirrors run_pathway): a seeded
  # current.msg would let CalculateHyperScore's error path segfault the process;
  # removing it forces the intended catchable R error, converted to a clean skip.
  if (exists("current.msg", envir = .GlobalEnv)) rm(list = "current.msg", envir = .GlobalEnv)
  mSet <- tryCatch(
    CalculateHyperScore(mSet),
    error = function(e)
      stop("ORA: too few mappable metabolites for enrichment (", conditionMessage(e), ")",
           call. = FALSE))
  if (!is.list(mSet)) {
    stop("ORA: too few mappable metabolites for enrichment", call. = FALSE)
  }

  # Expose the dir where CalculateHyperScore wrote msea_ora_result.csv (name-keyed).
  mSet$imic_workdir <- workdir
  mSet
}
