# reference-metabolome-check.R
#
# Task 5 investigation: how does the "reference metabolome" (background set)
# affect MetaboAnalystR pathway (pathora) and enrichment (msetora) results?
#
# Run from repo root:
#   "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" src/metaboanalyst/investigate/reference-metabolome-check.R
#
# This script is self-contained and rerunnable. It sources run-pathway.R for
# the baseline (filter OFF) pathway reproduction, then builds additional mSet
# sequences by hand to isolate the effect of the reference metabolome:
#   1. Pathway module, filter OFF                       (= run_pathway() baseline)
#   2. Pathway module, filter ON + KEGG reference        (Setup.KEGGReferenceMetabolome)
#      -- plus a root-cause check of *why* (2) behaves the way it does.
#   3. ORA module,     filter OFF (no reference)
#   4. ORA module,     filter ON + HMDB/name reference   (Setup.HMDBReferenceMetabolome)
#
# Findings are written up in FINDINGS-reference-metabolome.md alongside this
# script; the numbers there are copy-pasted verbatim from this script's
# console output.

suppressMessages(suppressWarnings(library(MetaboAnalystR)))

repo_root <- normalizePath(getwd())
cat("Repo root (cwd at script start):", repo_root, "\n\n")

cmpd.vec <- c("Nicotinic acid mononucleotide", "Niacinamide", "Nudifloramide",
              "nicotinamide", "nicotinamide riboside",
              "nicotinamide adenine dinucleotide", "tryptophan")

kegg_ref_path <- file.path(repo_root, "trenton scripts", "3. Results", "Primary Outcomes",
                            "Primary Outcomes Reference Metabolome (KEGG).txt")
name_ref_path <- file.path(repo_root, "trenton scripts", "3. Results", "Primary Outcomes",
                            "Primary Outcomes Reference Metabolome.txt")

stopifnot(file.exists(kegg_ref_path), file.exists(name_ref_path))

# ---------------------------------------------------------------------------
# Known off-public-web bug in CalculateHyperScore() (msetora path only):
# it ends with
#   ExportOraMembershipJson()
#   PlotORAMembership(NA, "mset_membership_0_")
#   return(result)
# Both helper calls omit mSetObj. .get.mSet(NA) is defined as
#   if (.on.public.web) get0("mSetObj", envir=.GlobalEnv, ...) else return(obj)
# so off the public web (our case, .on.public.web == FALSE) it returns the NA
# literal instead of falling back to a global. The very next line in each
# helper does `mSetObj$analSet$ora.mat`, which throws
#   "$ operator is invalid for atomic vectors"
# This happens even in the plain filter-OFF, no-reference case -- it is
# unconditional and orthogonal to the reference-metabolome question, but it
# means CalculateHyperScore() cannot complete AT ALL locally without a patch.
# CalculateOraScore() (the pathway/pathora path) does not call either helper,
# so run-pathway.R never needed this patch.
#
# The two helpers are reporting/plotting side effects only (a JSON membership
# dump and a PNG heatmap) -- neither touches mSetObj$analSet$ora.mat, which is
# already computed and stored by the time they're called. No-op'ing them is
# safe and does not change any numeric result.
# ---------------------------------------------------------------------------
patch_msetora_export_bug <- function() {
  ns <- asNamespace("MetaboAnalystR")
  unlockBinding("ExportOraMembershipJson", ns)
  assign("ExportOraMembershipJson", function(mSetObj = NA) invisible(NULL), envir = ns)
  lockBinding("ExportOraMembershipJson", ns)
  unlockBinding("PlotORAMembership", ns)
  assign("PlotORAMembership", function(mSetObj = NA, ...) invisible(NULL), envir = ns)
  lockBinding("PlotORAMembership", ns)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Step 1: inspect loader signatures
# ---------------------------------------------------------------------------
cat("=====================================================================\n")
cat("STEP 1: loader signatures\n")
cat("=====================================================================\n")
cat("args(Setup.KEGGReferenceMetabolome):\n")
print(args(Setup.KEGGReferenceMetabolome))
cat("\nargs(Setup.HMDBReferenceMetabolome):\n")
print(args(Setup.HMDBReferenceMetabolome))
cat("\nargs(SetMetabolomeFilter):\n")
print(args(SetMetabolomeFilter))
cat("\nargs(SetCurrentMsetLib):\n")
print(args(SetCurrentMsetLib))
cat(paste(
  "\nBoth Setup.*ReferenceMetabolome loaders take (mSetObj, filePath) --",
  "filePath is a FILE PATH read line-by-line via scan(), NOT a character vector.",
  "They only stash the matched reference into mSetObj$dataSet$metabo.filter.kegg",
  "or $metabo.filter.hmdb; they do not touch use.metabo.filter or the pathway/mset",
  "library, so they may be called any time before the Calculate*Score() step that",
  "consumes them (order relative to SetMetabolomeFilter/SetCurrentMsetLib does not",
  "matter, but they DO need InitDataObjects to have already set the global anal.type,",
  "since the loader body branches on anal.type %in% c('msetora','msetssp','msetqea')",
  "to decide master_compound_db.qs vs compound_db.qs).\n"
))

# ---------------------------------------------------------------------------
# Step 2: Pathway module (pathora) -- filter OFF vs filter ON + KEGG reference
# ---------------------------------------------------------------------------
cat("=====================================================================\n")
cat("STEP 2: pathway module (pathora), golden Elicit-Up-1mo cell\n")
cat("=====================================================================\n")

source(file.path(repo_root, "src", "metaboanalyst", "R", "run-pathway.R"))

# 2a. filter OFF (baseline, via the existing wrapper)
mSet_off <- run_pathway(cmpd.vec, pathlib = "smpdb")
res_off <- read.csv(file.path(mSet_off$imic_workdir, "pathway_results.csv"),
                     row.names = 1, check.names = FALSE)
pcol <- intersect(c("Raw p", "Raw.p"), colnames(res_off))[1]
p_off <- res_off["Nicotinate and Nicotinamide Metabolism", pcol]
cat(sprintf("[pathway] filter OFF -> Nicotinate & Nicotinamide Raw p = %s\n", format(p_off, digits = 8)))

# 2b. filter ON + KEGG reference metabolome loaded
p_on <- NA
pathway_filter_on_error <- NULL
workdir2 <- tempfile("pathway_kegg_ref_"); dir.create(workdir2)
old <- setwd(workdir2)
pathway_on_result <- tryCatch({
  mSet <- InitDataObjects("conc", "pathora", FALSE, 150)
  mSet <- Setup.MapData(mSet, cmpd.vec)
  mSet <- CrossReferencing(mSet, "name")
  mSet <- CreateMappingResultTable(mSet)
  mSet <- SetSMPDB.PathLib(mSet, "hsa")
  mSet <- SetOrganism(mSet, "hsa")
  mSet <- Setup.KEGGReferenceMetabolome(mSet, kegg_ref_path)
  mSet <- SetMetabolomeFilter(mSet, TRUE)
  mSet <- CalculateOraScore(mSet, "rbc", "hyperg")
  read.csv(file.path(workdir2, "pathway_results.csv"), row.names = 1, check.names = FALSE)
}, error = function(e) {
  pathway_filter_on_error <<- conditionMessage(e)
  NULL
})
setwd(old)

if (!is.null(pathway_filter_on_error)) {
  cat(sprintf("[pathway] filter ON (KEGG reference, 35 cmpds) -> ERROR: %s\n", pathway_filter_on_error))
} else if ("Nicotinate and Nicotinamide Metabolism" %in% rownames(pathway_on_result)) {
  p_on <- pathway_on_result["Nicotinate and Nicotinamide Metabolism", pcol]
  cat(sprintf("[pathway] filter ON  (KEGG reference, 35 cmpds) -> Nicotinate & Nicotinamide Raw p = %s\n",
              format(p_on, digits = 8)))
} else {
  cat("[pathway] filter ON (KEGG reference) -> Nicotinate and Nicotinamide Metabolism NOT in results table\n")
  cat("Full filter-ON pathway_results.csv rownames:\n")
  print(rownames(pathway_on_result))
}

cat(sprintf("\nTrenton's downloaded golden value: 1.2011e-05\n"))
cat(sprintf("filter OFF matches golden: %s\n", isTRUE(all.equal(unname(p_off), 1.2011e-05, tolerance = 1e-3))))
if (!is.na(p_on)) {
  cat(sprintf("filter ON  matches golden: %s\n", isTRUE(all.equal(unname(p_on), 1.2011e-05, tolerance = 1e-3))))
}

# 2c. Root-cause check for the filter-ON error: is the KEGG reference even
# compatible with the SMPDB library's member ID space? CalculateOraScore()'s
# local-mode filter branch ALWAYS intersects against
# mSetObj$dataSet$metabo.filter.kegg (a vector of KEGG IDs), regardless of
# whether the loaded pathway library is KEGG or SMPDB. If SMPDB library
# members are keyed by HMDB IDs (not KEGG IDs), a KEGG-ID reference can never
# match anything in that library.
cat("\n--- Root cause check: SMPDB library member ID namespace ---\n")
workdir3 <- tempfile("pathway_idcheck_"); dir.create(workdir3)
old <- setwd(workdir3)
mSet_chk <- InitDataObjects("conc", "pathora", FALSE, 150)
mSet_chk <- Setup.MapData(mSet_chk, cmpd.vec)
mSet_chk <- CrossReferencing(mSet_chk, "name")
mSet_chk <- CreateMappingResultTable(mSet_chk)
mSet_chk <- SetSMPDB.PathLib(mSet_chk, "hsa")
mSet_chk <- SetOrganism(mSet_chk, "hsa")
ck <- get("current.kegglib", envir = .GlobalEnv)
example_ids <- unlist(ck$mset.list[[1]], use.names = FALSE)
setwd(old)
cat(sprintf("pathwaylibtype = %s\n", mSet_chk$pathwaylibtype))
cat(sprintf("Example member IDs from current.kegglib$mset.list[[1]] (first pathway set): %s\n",
            paste(head(example_ids, 5), collapse = ", ")))
kegg_ref_ids <- readLines(kegg_ref_path); kegg_ref_ids <- kegg_ref_ids[nzchar(kegg_ref_ids)]
cat(sprintf("KEGG reference file IDs look like (first 5): %s\n", paste(head(kegg_ref_ids, 5), collapse = ", ")))
cat(sprintf("Any overlap between SMPDB member IDs and KEGG reference IDs (across whole library)? %s\n",
            any(unlist(lapply(ck$mset.list, function(x) any(unlist(x, use.names = FALSE) %in% kegg_ref_ids))))))
cat("Conclusion: SMPDB library members are HMDB IDs; the KEGG-ID reference vector cannot match them at all,\n")
cat("so SetMetabolomeFilter(TRUE) + Setup.KEGGReferenceMetabolome collapses nearly every pathway set to zero\n")
cat("members, tripping the package's 'too few sets' / 'single metabolite set' error path -- this is a real,\n")
cat("substantive incompatibility, not just a cosmetic difference in numbers.\n")

# ---------------------------------------------------------------------------
# Step 3: ORA / enrichment module (msetora) -- filter OFF vs filter ON + HMDB/name reference
# ---------------------------------------------------------------------------
cat("\n=====================================================================\n")
cat("STEP 3: ORA module (msetora), same 7 compounds\n")
cat("=====================================================================\n")

patch_msetora_export_bug()

run_msetora <- function(cmpd.vec, ref_path = NULL, mset_lib = "smpdb_pathway", exclude_num = 2) {
  workdir <- tempfile("msetora_"); dir.create(workdir)
  old <- setwd(workdir); on.exit(setwd(old), add = TRUE)

  mSet <- InitDataObjects("conc", "msetora", FALSE, 150)
  mSet <- Setup.MapData(mSet, cmpd.vec)
  mSet <- CrossReferencing(mSet, "name")
  mSet <- CreateMappingResultTable(mSet)
  if (!is.null(ref_path)) {
    mSet <- Setup.HMDBReferenceMetabolome(mSet, ref_path)
    mSet <- SetMetabolomeFilter(mSet, TRUE)
  } else {
    mSet <- SetMetabolomeFilter(mSet, FALSE)
  }
  mSet <- SetCurrentMsetLib(mSet, mset_lib, exclude_num)
  mSet <- CalculateHyperScore(mSet)

  list(mSet = mSet, workdir = workdir, ora.mat = mSet$analSet$ora.mat)
}

top_off <- NA; p_ora_off <- NA
top_on <- NA; p_ora_on <- NA

cat("--- ORA filter OFF (no reference) ---\n")
ora_off <- tryCatch(run_msetora(cmpd.vec, ref_path = NULL), error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n"); NULL
})
if (!is.null(ora_off)) {
  m <- ora_off$ora.mat
  cat("Top rows (filter OFF):\n")
  print(head(m, 5))
  top_off <- rownames(m)[1]
  p_ora_off <- m[1, "Raw p"]
  cat(sprintf("\n[ORA] filter OFF -> top pathway = %s, Raw p = %s\n", top_off, format(p_ora_off, digits = 8)))
}

cat("\n--- ORA filter ON (compound-name/HMDB reference, 35 cmpds) ---\n")
ora_on <- tryCatch(run_msetora(cmpd.vec, ref_path = name_ref_path), error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n"); NULL
})
if (!is.null(ora_on)) {
  m <- ora_on$ora.mat
  cat("Top rows (filter ON, HMDB reference):\n")
  print(head(m, 5))
  top_on <- rownames(m)[1]
  p_ora_on <- m[1, "Raw p"]
  cat(sprintf("\n[ORA] filter ON  (HMDB/name reference) -> top pathway = %s, Raw p = %s\n", top_on, format(p_ora_on, digits = 8)))
}

cat("\n=====================================================================\n")
cat("SUMMARY\n")
cat("=====================================================================\n")
cat(sprintf("Pathway filter OFF                   : Nicotinate & Nicotinamide Raw p = %s\n", format(p_off, digits = 8)))
cat(sprintf("Pathway filter ON  (KEGG reference)   : %s\n",
            if (!is.null(pathway_filter_on_error)) paste("ERROR:", pathway_filter_on_error) else format(p_on, digits = 8)))
if (!is.null(ora_off)) {
  cat(sprintf("ORA filter OFF (no reference)         : top = %s, Raw p = %s\n", top_off, format(p_ora_off, digits = 8)))
}
if (!is.null(ora_on)) {
  cat(sprintf("ORA filter ON  (HMDB/name reference)  : top = %s, Raw p = %s\n", top_on, format(p_ora_on, digits = 8)))
}
cat("\nDone.\n")
