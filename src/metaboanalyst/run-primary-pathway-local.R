# run-primary-pathway-local.R -- web-independent Figure 3B (primary KEGG pathway-impact).
#
# Figure 3B plots pathway Impact (x) vs -log10(Raw p) (y) per study x timepoint.
# It was pinned to Trenton's DOWNLOADED MetaboAnalyst pathway_results.csv exports
# (results/metaboanalyst/primary_pathway_trenton/) because an earlier local run used
# the SMPDB library (larger sets) which COMPRESSED the impact axis.
#
# This runner reproduces the submitted exports LOCALLY by running the SAME KEGG
# library the web tool uses (SetKEGG.PathLib "metpa"; src/metaboanalyst/R/run-pathway.R
# with pathlib="kegg"). Verified byte-identical to Trenton's export on the Elicit
# 1-month cell: Nicotinate and nicotinamide metabolism Total=15, Hits=5,
# Raw p=3.5461e-11, Impact=0.61974 -- so 3B no longer needs any web download.
#
# Cells = Trenton's Fig 3B construction: one non-directional cell per study x
# timepoint, query = combined-arm sigFDR==1 primary metabolites (pos + neg together),
# via build_cells_combined_sigfdr() -- identical to run-primary-pathway-compare.R.
#
# Output: results/metaboanalyst/primary_pathway_local/primary_pathway_all_cells.csv
#         (study, tp, pathway, total, expected, hits, impact, raw_p, fdr)
# Run from repo root:  Rscript src/metaboanalyst/run-primary-pathway-local.R
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
