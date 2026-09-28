# =============================================================================
# run-pathway.R
#
# Helpers for MetaboAnalystR's pathway module (pathora). run_pathway() runs one
# query list against the KEGG or SMPDB human library with a hypergeometric test
# and relative-betweenness topology, metabolome filter off (no reference
# metabolome), following the call sequence of the saved metaboanalyst.ca pathway
# analysis, which it reproduces exactly. harvest_pathway() turns a KEGG run into
# a pathway table and a compound -> pathway hit table. Sourced by
# run-primary-pathway-local.R (Fig 3B, Table S2) and run-compartment-pathway.R
# (Fig 6D, Table S11).
#
# Inputs : src/metaboanalyst/reference/kegg_hsa_pathway_names.csv (KEGG id -> pathway name)
# Outputs: none (returns an mSet / tables); MetaboAnalystR downloads its libraries
#          from metaboanalyst.ca and they are cached per R session
# =============================================================================

run_pathway <- function(cmpd.vec,
                        pathlib     = c("smpdb", "kegg"),
                        topology    = "rbc",
                        enrich      = "hyperg",
                        kegg.version = "current",
                        dpi         = 150) {
  pathlib <- match.arg(pathlib)
  stopifnot(length(cmpd.vec) >= 1)

  # MetaboAnalystR writes result CSVs to the working directory, so give each run
  # its own throwaway dir and restore the caller's wd on exit.
  workdir <- tempfile("pathway_"); dir.create(workdir)
  original_wd <- setwd(workdir); on.exit(setwd(original_wd), add = TRUE)

  # Library cache, as in run_ora(): MetaboAnalystR re-downloads its libraries into
  # every fresh workdir, and a failed download silently drops the cell. One cache per
  # pathway library, because KEGG and SMPDB both name theirs "hsa.qs";
  # current.kegglib.qs is per-run state and is never cached.
  lib_files <- c("compound_db.qs", "syn_nms.qs", "hsa.qs")
  lib_cache <- file.path(tempdir(), "metaboanalyst_libs", pathlib)
  dir.create(lib_cache, recursive = TRUE, showWarnings = FALSE)
  cached <- intersect(list.files(lib_cache), lib_files)
  file.copy(file.path(lib_cache, cached), workdir, copy.date = TRUE)

  # dpi arg is mandatory in MetaboAnalystR 4.3.0 (self-referential default bug).
  mSet <- InitDataObjects("conc", "pathora", FALSE, dpi)
  mSet <- Setup.MapData(mSet, cmpd.vec)
  mSet <- CrossReferencing(mSet, "name")
  mSet <- CreateMappingResultTable(mSet)
  # SetKEGG.PathLib in MetaboAnalystR 4.3.0 requires a lib.version arg ("current"
  # = the offline KEGG metpa library the package downloads once and caches, same
  # one-time DB fetch as the compound DB; no per-analysis web-tool call).
  if (pathlib == "smpdb") mSet <- SetSMPDB.PathLib(mSet, "hsa")
  else                    mSet <- SetKEGG.PathLib(mSet, "hsa", kegg.version)
  mSet <- SetOrganism(mSet, "hsa")
  mSet <- SetMetabolomeFilter(mSet, FALSE)   # filter off = reproduces the paper's numbers

  # Robust too-few-metabolites handling. When fewer than 3 compounds map into the
  # pathway library, MetaboAnalystR's CalculateOraScore takes an error path via
  # AddErrMsg(). If `current.msg` happens to be seeded in the global env (e.g. by
  # an earlier run_ora() in the same session), that path proceeds into C code that
  # segfaults the whole process. Removing the seed forces the intended catchable R
  # error instead, which we convert into a clean, informative skip reason.
  if (exists("current.msg", envir = .GlobalEnv)) rm(list = "current.msg", envir = .GlobalEnv)
  mSet <- tryCatch(
    CalculateOraScore(mSet, topology, enrich),
    error = function(e)
      stop("pathway analysis: too few mappable metabolites for enrichment (", conditionMessage(e), ")",
           call. = FALSE)
  )
  if (!is.list(mSet)) {
    stop("pathway analysis: too few mappable metabolites for enrichment", call. = FALSE)
  }

  new_libs <- setdiff(intersect(list.files(workdir), lib_files), cached)
  file.copy(file.path(workdir, new_libs), lib_cache, copy.date = TRUE)

  # Expose the dir where CalculateOraScore wrote pathway_results.csv (name-keyed).
  mSet$imic_workdir <- workdir
  mSet
}

# harvest_pathway() -- turn a run_pathway() KEGG mSet into two tidy tables, in
# place of the web tool's downloaded results and a hand-coded compound -> pathway map:
#   $pathways  <- the "pathway_results.csv" content (KEGG topology + ORA stats),
#                 with readable pathway names (the metpa hsa-id -> name map, bundled
#                 at src/metaboanalyst/reference/kegg_hsa_pathway_names.csv so the
#                 labels match the web tool exactly, e.g. "Citrate cycle (TCA cycle)").
#   $metabolite_pathways <- compound -> pathway rows from the ORA hit membership:
#                 for each pathway, ora.hits gives the hit KEGG compound ids, mapped
#                 back to the original input names via the name map.table
#                 (Query -> KEGG). One row per (input compound, pathway).
# KEGG's ora.mat is keyed by hsa id (this MetaboAnalystR build does not write readable
# names into the CSV), hence the bundled id->name lookup.
KEGG_NAME_MAP <- "src/metaboanalyst/reference/kegg_hsa_pathway_names.csv"

harvest_pathway <- function(mSet, name_map_file = KEGG_NAME_MAP) {
  id2name <- { m <- utils::read.csv(name_map_file, stringsAsFactors = FALSE)
               stats::setNames(m$pathway, m$kegg_id) }
  om <- as.data.frame(mSet$analSet$ora.mat, check.names = FALSE)
  om$kegg_id <- rownames(om)
  om$pathway <- ifelse(is.na(id2name[om$kegg_id]), om$kegg_id, id2name[om$kegg_id])
  pathways <- data.frame(
    pathway = om$pathway, kegg_id = om$kegg_id,
    total = om$Total, expected = om$Expected, hits = om$Hits,
    raw_p = om[["Raw p"]], fdr = om$FDR, impact = om$Impact,
    stringsAsFactors = FALSE)
  pathways <- pathways[order(pathways$raw_p), ]

  # metabolite -> pathway from ORA hit membership
  mt        <- as.data.frame(mSet$dataSet$map.table, stringsAsFactors = FALSE)
  kegg2query <- stats::setNames(mt$Query, mt$KEGG)
  oh <- mSet$analSet$ora.hits; oh <- oh[lengths(oh) > 0]
  mp <- do.call(rbind, lapply(names(oh), function(pid) {
    ids <- unname(oh[[pid]]); qn <- kegg2query[ids]; qn <- qn[!is.na(qn)]
    if (!length(qn)) return(NULL)
    data.frame(tracking_name = unname(qn),
               pathway = ifelse(is.na(id2name[pid]), pid, id2name[[pid]]),
               kegg_id = pid,
               raw_p = mSet$analSet$ora.mat[pid, "Raw p"],
               fdr   = mSet$analSet$ora.mat[pid, "FDR"],
               impact = mSet$analSet$ora.mat[pid, "Impact"],
               stringsAsFactors = FALSE)
  }))
  if (!is.null(mp)) mp <- mp[order(mp$raw_p), ]
  list(pathways = pathways, metabolite_pathways = mp)
}
