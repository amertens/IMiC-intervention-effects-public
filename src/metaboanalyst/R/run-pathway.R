# run-pathway.R — thin wrapper around MetaboAnalystR's pathora engine.
# Reproduces the sequence Task 2 verified against Trenton's golden pathway cell.
# Runs filter OFF (no reference metabolome), matching the paper.

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
  # SEGFAULTS the whole process. Removing the seed forces the intended catchable R
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

  # Expose the dir where CalculateOraScore wrote pathway_results.csv (name-keyed).
  mSet$imic_workdir <- workdir
  mSet
}

# harvest_pathway() -- turn a run_pathway() KEGG mSet into two tidy, web-independent
# tables that replace Trenton's two web/manual inputs in Compartment Tracking.Rmd:
#   $pathways  <- the downloaded "pathway_results.csv" (KEGG topology + ORA stats),
#                 with readable pathway names (the metpa hsa-id -> name map, bundled
#                 at src/metaboanalyst/reference/kegg_hsa_pathway_names.csv so the
#                 labels match the web tool exactly, e.g. "Citrate cycle (TCA cycle)").
#   $metabolite_pathways <- the hardcoded metabolite->pathway tribble, reconstructed
#                 from the ORA hit membership: for each pathway, ora.hits gives the hit
#                 KEGG compound ids, mapped back to the ORIGINAL input names via the
#                 name map.table (Query -> KEGG). One row per (input compound, pathway).
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
