# =============================================================================
# run-primary-pathway-local.R
#
# KEGG pathway analysis of the primary milk outcomes for Fig 3B and Table S2.
# One non-directional cell per study x time point x contrast: the query is the
# combined-arm primary metabolites with FDR < 0.05 (up- and down-regulated
# together, build_cells_combined_sigfdr()), run through MetaboAnalystR's pathway
# module with the KEGG (metpa) library, hypergeometric test and
# relative-betweenness topology, no reference-metabolome filter
# (R/run-pathway.R). Fig 3B plots pathway impact against -log10(raw P).
# The KEGG library matches the web tool's; on the ELICIT 1-month cell this run
# reproduced the saved web export (nicotinate and nicotinamide metabolism:
# total 15, hits 5, raw P 3.5461e-11, impact 0.61974).
#
# Inputs : [needs restricted data] results/combined_intervention_effects_results_combined_arms.RDS
#          src/metaboanalyst/reference/kegg_hsa_pathway_names.csv  (via R/run-pathway.R)
# Outputs: results/metaboanalyst/primary_pathway_local/primary_pathway_all_cells.csv   (Fig 3B, Table S2)
#          results/metaboanalyst/primary_pathway_local/primary_pathway_hits_all_cells.csv (compound -> pathway hits)
#
# MetaboAnalystR downloads its compound and KEGG libraries from metaboanalyst.ca
# at run time.
# Run from repo root: Rscript src/metaboanalyst/run-primary-pathway-local.R
# =============================================================================
suppressMessages({ library(dplyr); library(stringr); library(readr); library(MetaboAnalystR) })
source("src/metaboanalyst/R/build-cells.R")
source("src/metaboanalyst/R/run-pathway.R")

COMBINED_RDS <- "results/combined_intervention_effects_results_combined_arms.RDS"
OUT_DIR      <- "results/metaboanalyst/primary_pathway_local"

# "Misame (14-21 days)" -> study "MISAME-III", tp "14-21 days" (matches fig study names)
.study_from <- function(s) dplyr::case_when(
  grepl("Elicit", s) ~ "ELICIT", grepl("Misame", s) ~ "MISAME-III",
  grepl("Vital",  s) ~ "Mumta-LW", TRUE ~ s)
.tp_from <- function(s) trimws(gsub("[()]", "", str_extract(s, "\\(([^)]+)\\)")))

if (sys.nframe() == 0) {
  combined <- readRDS(COMBINED_RDS)
  cells <- build_cells_combined_sigfdr(combined, "primary")
  message("primary pathway cells (study x timepoint): ", length(cells))

  rows <- list()
  hitrows <- list()
  for (cell in cells) {
    mSet <- tryCatch(run_pathway(cell$query, pathlib = "kegg"),
                     error = function(e) { message("  skip ", cell$study, ": ", conditionMessage(e)); NULL })
    if (is.null(mSet)) next
    h <- harvest_pathway(mSet)
    st <- .study_from(cell$study); tpv <- .tp_from(cell$study)
    rows[[length(rows) + 1]] <- h$pathways %>%
      mutate(study = st, tp = tpv, .before = 1)
    if (!is.null(h$metabolite_pathways) && nrow(h$metabolite_pathways)) {
      hitrows[[length(hitrows) + 1]] <- h$metabolite_pathways %>%
        mutate(study = st, tp = tpv, .before = 1) %>%
        select(study, tp, pathway, kegg_id, tracking_name)
    }
    message("  ", cell$study, ": ", nrow(h$pathways), " pathways, ",
            sum(h$pathways$raw_p < 0.05), " at p<0.05")
  }
  tab <- bind_rows(rows) %>%
    select(study, tp, pathway, kegg_id, total, expected, hits, impact, raw_p, fdr)
  hits <- bind_rows(hitrows)

  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  write_csv(tab, file.path(OUT_DIR, "primary_pathway_all_cells.csv"))
  write_csv(hits, file.path(OUT_DIR, "primary_pathway_hits_all_cells.csv"))
  message("wrote ", nrow(tab), " rows -> ", OUT_DIR,
          " | impact range: ", paste(round(range(tab$impact, na.rm = TRUE), 3), collapse = "-"))
  message("wrote ", nrow(hits), " hit rows -> primary_pathway_hits_all_cells.csv")
}
