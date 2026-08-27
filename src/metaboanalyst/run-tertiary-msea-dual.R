# run-tertiary-msea-dual.R -- fully-automated Fig 5B reproduction that mirrors the
# MetaboAnalyst web workflow's TWO runs and unions them:
#   * metabolite pass: amino-acid/energy compounds vs the metabolite reference,
#     standard compound_db name matching (run_ora lipid = FALSE);
#   * lipid pass: Quant 500 lipids reformatted to LIPID MAPS names (see
#     R/lipid-name-map.R) vs the lipid reference, lipid_compound_db matching
#     (run_ora lipid = TRUE).
# Each class runs in a SEPARATE Rscript subprocess so MetaboAnalystR's global
# compound-DB state (compound_db vs lipid_compound_db) does not leak between passes
# -- mixing them in one session silently breaks name matching.
#
# Usage (from repo root):
#   Rscript src/metaboanalyst/run-tertiary-msea-dual.R              # driver: both passes + union
#   Rscript src/metaboanalyst/run-tertiary-msea-dual.R metabolite  # single pass (internal)
#   Rscript src/metaboanalyst/run-tertiary-msea-dual.R lipid
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
  # Shared background = Andrew's saved reference metabolome (metabolites+lipids,
  # ~1,268 names). The submitted 5B used ONE reference for BOTH feature-type runs;
  # the class split is on the QUERY only. Convert the reference lipids to LIPID MAPS
  # for the lipid pass so they match the lipid DB.
  ref_qer <- unique(trimws(readLines("src/metaboanalyst/reference/refMetabolomeForQER.csv", warn = FALSE)))
  ref_qer <- ref_qer[ref_qer != ""]
  ref <- if (is_lip) unique(to_lipidmaps(ref_qer)) else ref_qer
  message(sprintf("[%s] reference (QER shared background): %d names", klass, length(ref)))

  rows <- list()
  # Union combined-arm and stratified-arm cells: the submitted from-tables mixes
  # combined ELICIT/Mumta contrasts with stratified MISAME arms. For the lipid pass
  # rewrite lipid label_f -> LIPID MAPS; restrict tertiary rows to this class
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
                       error = function(e) NULL)
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
    apply_pathway_size_floor(min_size = 1) %>% arrange(raw_p)  # submitted from-tables kept total=1 lipid pathways
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
            file.path(OUT, "tertiary_msea_dual_perContrastFDRsig_tableS4.csv"))

  # Compare recovery to the submitted from-tables panel.
  sub <- tryCatch(read_csv(file.path(OUT, "tertiary_msea_fromTables.csv"), show_col_types = FALSE),
                  error = function(e) NULL)
  dual_sig <- unique(dual$pathway[dual$fdr_native < 0.05 & !is.na(dual$fdr_native)])
  message(sprintf("\n=== DUAL union: %d FDR-sig rows | %d distinct pathways (metab %d + lipid %d) ===",
    sum(dual$fdr_native < 0.05, na.rm = TRUE), length(dual_sig),
    dplyr::n_distinct(dual$pathway[dual$klass=="metabolite" & dual$fdr_native<0.05 & !is.na(dual$fdr_native)]),
    dplyr::n_distinct(dual$pathway[dual$klass=="lipid"      & dual$fdr_native<0.05 & !is.na(dual$fdr_native)])))
  if (!is.null(sub)) {
    sub_sig <- unique(sub$pathway[sub$fdr_native < 0.05 & !is.na(sub$fdr_native)])
    message(sprintf("submitted from-tables distinct FDR-sig pathways: %d", length(sub_sig)))
    message("recovered (in both): ", paste(sort(intersect(dual_sig, sub_sig)), collapse = " | "))
    message("MISSED (submitted, not dual): ", paste(sort(setdiff(sub_sig, dual_sig)), collapse = " | "))
    message(sprintf("RECOVERY: %d / %d submitted pathways (%.0f%%)",
      length(intersect(dual_sig, sub_sig)), length(sub_sig),
      100 * length(intersect(dual_sig, sub_sig)) / max(1, length(sub_sig))))
  }
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) && args[1] %in% c("metabolite", "lipid")) run_one_pass(args[1]) else driver()
}
