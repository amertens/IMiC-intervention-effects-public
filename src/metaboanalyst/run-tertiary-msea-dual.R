# =============================================================================
# run-tertiary-msea-dual.R
#
# Over-representation analysis (ORA) of the tertiary (targeted metabolome) milk
# outcomes for Fig 5B and Table S3. The MetaboAnalyst workflow runs this as two
# feature types and pools them:
#   - metabolite pass: non-lipid compounds, standard compound_db name matching
#     (run_ora(lipid = FALSE));
#   - lipid pass: Quant 500 lipids renamed to LIPID MAPS abbreviations
#     (R/lipid-name-map.R), lipid_compound_db name matching (run_ora(lipid = TRUE)).
# Cells are directional (up/down, P < 0.05) per study x time point x contrast,
# from both the combined-arm and the stratified-arm results; the background is
# the shared 1,268-name reference metabolome for both passes.
# Each pass runs in its own Rscript subprocess because MetaboAnalystR keeps the
# compound database (compound_db vs lipid_compound_db) in global state, and
# mixing the two in one session breaks name matching.
#
# Inputs : [needs restricted data]
#          results/combined_intervention_effects_results_combined_arms.RDS
#          results/combined_intervention_effects_results_stratified_arms.RDS
#          src/metaboanalyst/reference/refMetabolomeForQER.csv  (reference metabolome)
# Outputs: results/metaboanalyst/tertiary_msea/tertiary_msea_dual.csv          (Fig 5B, Table S3)
#          results/metaboanalyst/tertiary_msea/tertiary_msea_dual_fdr_sig.csv  (fdr_native < 0.05 subset)
#          results/metaboanalyst/tertiary_msea/tertiary_msea_dual_{metabolite,lipid}.csv (per-pass tables)
#
# MetaboAnalystR downloads its compound and metabolite-set libraries from
# metaboanalyst.ca at run time.
# Usage (from repo root):
#   Rscript src/metaboanalyst/run-tertiary-msea-dual.R              # driver: both passes + union
#   Rscript src/metaboanalyst/run-tertiary-msea-dual.R metabolite  # single pass (called by the driver)
#   Rscript src/metaboanalyst/run-tertiary-msea-dual.R lipid
# =============================================================================
suppressMessages({ library(dplyr); library(readr); library(stringr) })
source("src/metaboanalyst/R/label-map.R")
source("src/metaboanalyst/R/build-cells.R")
source("src/metaboanalyst/R/run-ora.R")
source("src/metaboanalyst/R/harvest.R")
source("src/metaboanalyst/R/build-reference.R")
source("src/metaboanalyst/R/lipid-name-map.R")

COMBINED_RDS   <- "results/combined_intervention_effects_results_combined_arms.RDS"
STRATIFIED_RDS <- "results/combined_intervention_effects_results_stratified_arms.RDS"
OUT <- "results/metaboanalyst/tertiary_msea"

.study_label <- function(st) {
  base <- str_replace(st, "\\s*\\(.*$", "")
  dplyr::recode(base, "Misame"="MISAME-III", "Vital"="Mumta-LW", "Elicit"="ELICIT", .default=base)
}
.timepoint <- function(st) str_trim(str_replace(str_replace(st, "^[^(]*\\(", ""), "\\)\\s*$", ""))

run_one_pass <- function(klass) {
  is_lip <- klass == "lipid"
  # One shared background (metabolites + lipids, ~1,268 names) for both feature
  # types; the class split applies to the query only. The reference lipids are
  # converted to LIPID MAPS names for the lipid pass so they match the lipid DB.
  ref_qer <- unique(trimws(readLines("src/metaboanalyst/reference/refMetabolomeForQER.csv", warn = FALSE)))
  ref_qer <- ref_qer[ref_qer != ""]
  ref <- if (is_lip) unique(to_lipidmaps(ref_qer)) else ref_qer
  message(sprintf("[%s] reference (QER shared background): %d names", klass, length(ref)))

  rows <- list()
  # Pool combined-arm and stratified-arm cells: the published table uses combined
  # ELICIT/Mumta-LW contrasts and stratified MISAME-III arms. For the lipid pass,
  # rewrite lipid label_f to LIPID MAPS; restrict tertiary rows to this class
  # (apply_label_map inside build_cells leaves LIPID MAPS names untouched).
  for (rds in c(COMBINED_RDS, STRATIFIED_RDS)) {
    if (!file.exists(rds)) next
    dat <- readRDS(rds) %>%
      mutate(label_f = if_else(is_lip & outcome_group == "tertiary" & category %in% LIPID_CATEGORIES,
                               to_lipidmaps(label_f), label_f)) %>%
      filter(!(outcome_group == "tertiary" &
               ((is_lip & !(category %in% LIPID_CATEGORIES)) |
                (!is_lip &  (category %in% LIPID_CATEGORIES)))))
    for (cell in build_cells(dat, "tertiary")) {
      mSet <- tryCatch(run_ora(cell$query, reference_names = ref, lipid = is_lip),
                       error = function(e) {
                         message("  skipped ", cell$study, " / ", cell$direction, ": ", conditionMessage(e))
                         NULL
                       })
      if (is.null(mSet)) next
      res <- tryCatch(harvest_results(mSet), error = function(e) NULL)
      if (is.null(res) || !nrow(res)) next
      res <- res %>% mutate(
        enrichment_ratio = ifelse(cell$direction == "down", -1, 1) * (hits / expected),
        study = .study_label(cell$study), timepoint = .timepoint(cell$study),
        contrast = cell$contrast, direction = cell$direction, klass = klass)
      rows[[length(rows) + 1]] <- res
    }
  }
  tab <- bind_rows(rows) %>%
    transmute(study, timepoint, contrast, direction, klass, pathway,
              total, expected, hits, raw_p, fdr_native = fdr, enrichment_ratio) %>%
    apply_pathway_size_floor(min_size = 1) %>% arrange(raw_p)  # no size floor: single-member lipid pathways are kept
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  write_csv(tab, file.path(OUT, paste0("tertiary_msea_dual_", klass, ".csv")))
  message(sprintf("[%s] wrote %d rows | %d FDR-sig | %d distinct sig pathways",
    klass, nrow(tab), sum(tab$fdr_native < 0.05, na.rm = TRUE),
    dplyr::n_distinct(tab$pathway[tab$fdr_native < 0.05 & !is.na(tab$fdr_native)])))
  invisible(tab)
}

driver <- function() {
  Rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  self <- "src/metaboanalyst/run-tertiary-msea-dual.R"
  for (k in c("metabolite", "lipid")) {
    message("=== subprocess pass: ", k, " ===")
    system2(Rscript, c(shQuote(self), k))
  }
  met <- read_csv(file.path(OUT, "tertiary_msea_dual_metabolite.csv"), show_col_types = FALSE)
  lip <- read_csv(file.path(OUT, "tertiary_msea_dual_lipid.csv"),      show_col_types = FALSE)
  dual <- bind_rows(met, lip) %>% arrange(raw_p)
  write_csv(dual, file.path(OUT, "tertiary_msea_dual.csv"))
  write_csv(filter(dual, fdr_native < 0.05),
            file.path(OUT, "tertiary_msea_dual_fdr_sig.csv"))

  dual_sig <- unique(dual$pathway[dual$fdr_native < 0.05 & !is.na(dual$fdr_native)])
  message(sprintf("\n=== DUAL union: %d FDR-sig rows | %d distinct pathways (metab %d + lipid %d) ===",
    sum(dual$fdr_native < 0.05, na.rm = TRUE), length(dual_sig),
    dplyr::n_distinct(dual$pathway[dual$klass=="metabolite" & dual$fdr_native<0.05 & !is.na(dual$fdr_native)]),
    dplyr::n_distinct(dual$pathway[dual$klass=="lipid"      & dual$fdr_native<0.05 & !is.na(dual$fdr_native)])))
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  # without the feature-level result files every subprocess pass fails, and the driver
  # would still exit 0; stop with the reason instead
  if (!any(file.exists(c(COMBINED_RDS, STRATIFIED_RDS))))
    stop("needs ", COMBINED_RDS, " / ", STRATIFIED_RDS,
         " (feature-level result files, available from the authors on request)", call. = FALSE)
  if (length(args) && args[1] %in% c("metabolite", "lipid")) run_one_pass(args[1]) else driver()
}
