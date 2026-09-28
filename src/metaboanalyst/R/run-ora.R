# =============================================================================
# run-ora.R
#
# Helper: run_ora() wraps MetaboAnalystR's over-representation module (msetora)
# for one query list: name matching (compound or lipid database), an optional
# reference metabolome as the background (metabolome filter on), then the
# hypergeometric test against a metabolite-set library (default SMPDB pathways).
# It includes the session patch for a MetaboAnalystR 4.3.0 crash in
# CalculateHyperScore() off the public web (patch_msetora_export_bug() below).
# Sourced by run-tertiary-msea-dual.R (Fig 5B, Table S3) and
# run-untargeted-msea.R (Fig 6A, Table S5).
#
# Inputs : none (takes character vectors)
# Outputs: none (returns an mSet); MetaboAnalystR downloads its libraries from
#          metaboanalyst.ca and they are cached per R session
# =============================================================================

# ---------------------------------------------------------------------------
# Known off-public-web bug in CalculateHyperScore() (msetora path only): it
# ends with ExportOraMembershipJson() and PlotORAMembership(NA, ...), both of
# which omit mSetObj. .get.mSet(NA) short-circuits to the NA literal off the
# public web (instead of falling back to a global mSetObj), so the very next
# line in each helper (`mSetObj$analSet$ora.mat`) throws "$ operator is
# invalid for atomic vectors". This happens unconditionally, even in the
# plain filter-off/no-reference case. Both helpers are reporting/plotting
# side effects only (a JSON membership dump and a PNG heatmap) that run after
# ora.mat is already computed and stored, so no-op'ing them does not change
# the numbers (ora.mat is identical with and without the patch).
#
# AddErrMsg() also reads an uninitialized `current.msg` global off the public
# web. It is not seeded before scoring: a seeded value lets the
# too-few-metabolites error path proceed into C code that segfaults the process
# (same failure handled in run_pathway). Removing current.msg before
# CalculateHyperScore makes that path raise a catchable R error, which is
# converted to a clean skip reason.
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

  # MetaboAnalystR downloads its libraries from metaboanalyst.ca into the working
  # directory and reuses them only if they are already there, so a fresh workdir per
  # call would re-download ~6 MB on every call, plus the 12.5 MB master_compound_db.qs
  # whenever a reference metabolome is set (slow, and a failed download drops the
  # cell). Seed each workdir from a per-session cache and save new downloads back
  # once the run succeeds. Whitelisted names only: current.msetlib.qs is per-run state.
  lib_files <- c("compound_db.qs", "syn_nms.qs", "lipid_compound_db.qs", "lipid_syn_nms.qs",
                 "master_compound_db.qs", paste0(mset_lib, ".qs"))
  lib_cache <- file.path(tempdir(), "metaboanalyst_libs", "ora")
  dir.create(lib_cache, recursive = TRUE, showWarnings = FALSE)
  cached <- intersect(list.files(lib_cache), lib_files)
  file.copy(file.path(lib_cache, cached), workdir, copy.date = TRUE)

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

  new_libs <- setdiff(intersect(list.files(workdir), lib_files), cached)
  file.copy(file.path(workdir, new_libs), lib_cache, copy.date = TRUE)

  # Expose the dir where CalculateHyperScore wrote msea_ora_result.csv (name-keyed).
  mSet$imic_workdir <- workdir
  mSet
}
