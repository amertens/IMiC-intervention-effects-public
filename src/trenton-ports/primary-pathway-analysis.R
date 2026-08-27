# =============================================================================
# primary-pathway-analysis.R
# Faithful R port of Trenton's "Primary Outcomes (Pathway Analysis).Rmd".
#
# Mirrors his notebook:
#   1. Import the COMBINED-arm intervention-effects RDS (this analysis is combined
#      arms only -- no stratified cells, matching his Rmd).
#   2. Generate the MetaboAnalyst compound label (his verbatim 11-entry synonym map).
#   3. Per (study x timepoint) cell: query = primary features with sigFDR == 1,
#      "Total *" roll-ups excluded, POSITIVE AND NEGATIVE estimates analysed TOGETHER
#      (his rule, line 76: "Positive and negative estimates are analyzed together").
#      Run Pathway Analysis (pathora: hypergeometric + relative-betweenness topology,
#      SMPDB / Homo sapiens).
#
# AUTOMATION: his manual metaboanalyst.ca "Pathway Analysis" upload/download is
# replaced by run_pathway(), which reproduces his ONE downloaded pathway cell
# (Elicit Up 1mo) EXACTLY -- 6/6 shared pathways, relative difference 0.
#
# ONE DOCUMENTED FIDELITY DECISION (not our error, and not a change to his result):
#   His written instructions upload a KEGG reference metabolome for the pathway
#   step, and his Rmd note says "the same reference metabolome is used". However,
#   his DOWNLOADED pathway_results.csv reproduces ONLY with the metabolome filter
#   OFF (no reference): filter ON with a KEGG-keyed reference is broken on the
#   SMPDB library (its set members are HMDB IDs -> zero overlap -> MetaboAnalyst's
#   "too few sets" error). So run_pathway() runs filter OFF, which matches his
#   actual downloaded numbers. (This is the open "3B reference" item; the port
#   reproduces his ground truth rather than his written-but-ineffective upload.)
#
# Run from repo root:  Rscript src/trenton-ports/primary-pathway-analysis.R
# =============================================================================
suppressMessages({ library(dplyr); library(readr); library(stringr); library(MetaboAnalystR) })
source("src/metaboanalyst/R/label-map.R")        # his verbatim synonym map
source("src/metaboanalyst/R/build-cells.R")      # build_cells_combined_sigfdr(): his sigFDR==1, pos+neg-together rule
source("src/metaboanalyst/R/run-pathway.R")      # run_pathway(): scripted Pathway Analysis (pathora, filter OFF)
source("src/metaboanalyst/R/harvest.R")          # harvest_results(): tidy the downloaded result (incl. Impact)

COMBINED_RDS <- "trenton scripts/1. Data/combined_intervention_effects_results_combined_arms.RDS"
OUT          <- "results/trenton-ports/primary_pathway"

.cohort <- function(studytime) {
  s <- str_replace(studytime, "\\s*\\(.*$", "")
  recode(s, "Elicit" = "ELICIT", "Misame" = "MISAME-III", "Vital" = "Mumta-LW",
         "Mumpta" = "Mumta-LW", .default = s)
}
.timepoint <- function(studytime) str_trim(str_replace(str_replace(studytime, "^[^(]*\\(", ""), "\\)\\s*$", ""))

if (sys.nframe() == 0) {
  combined <- readRDS(COMBINED_RDS)

  # His "Analysis": one non-directional sigFDR==1 cell per study x timepoint
  # (pos + neg together). build_cells_combined_sigfdr() encodes exactly this.
  cells <- build_cells_combined_sigfdr(combined, "primary")
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  for (cell in cells) {
    mSet <- tryCatch(run_pathway(cell$query), error = function(e) { message("  skip: ", cell$study, " (", conditionMessage(e), ")"); NULL })
    if (is.null(mSet)) next
    res <- harvest_results(mSet)   # pathway, total, hits, expected, raw_p, fdr, impact
    if (!nrow(res)) next
    res <- res %>% mutate(study = .cohort(cell$study), timepoint = .timepoint(cell$study),
                          contrast = cell$contrast, .before = 1)
    slug <- gsub("[^A-Za-z0-9]+", "_", cell$study)
    celldir <- file.path(OUT, slug); dir.create(celldir, showWarnings = FALSE)
    write_csv(res, file.path(celldir, "pathway_results.csv"))
    rows[[length(rows) + 1]] <- res
  }

  if (!length(rows)) stop("no pathway cells produced results -- check MetaboAnalystR is installed and cells are non-empty")
  all <- bind_rows(rows) %>% arrange(study, timepoint, raw_p)
  write_csv(all, file.path(OUT, "primary_pathway_all_cells.csv"))
  message(sprintf("done: %d pathway rows across %d cells; %d FDR-significant. -> %s",
                  nrow(all), nrow(distinct(all, study, timepoint)),
                  sum(all$fdr < 0.05, na.rm = TRUE), OUT))
}
