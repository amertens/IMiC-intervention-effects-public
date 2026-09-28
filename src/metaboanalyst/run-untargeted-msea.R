# =============================================================================
# run-untargeted-msea.R
#
# Over-representation analysis (ORA) of the untargeted milk metabolome for Fig 6A
# and Table S5. For each study x time point x contrast x direction cell, the query
# is the annotated names of features with unadjusted P < 0.05, tested with
# MetaboAnalystR against the SMPDB pathway metabolite sets, using the 615-name
# reference metabolome as the background. Pathways with fewer than three members
# are dropped from the reported table; fdr_native is MetaboAnalyst's per-cell BH FDR.
#
# Inputs : results/milk_nominal_putative_annotation.csv  (src/2 analysis/47-annotate-milk-features.R)
#          src/metaboanalyst/reference/reference_metabolome_1593_matched_names.txt
#            (615 names; built by build-reference-metabolome.R)
# Outputs: results/metaboanalyst/untargeted_msea/untargeted_msea_combined.csv          (Fig 6A)
#          results/metaboanalyst/untargeted_msea/untargeted_msea_combined_fdr_sig.csv  (Table S5; fdr_native < 0.05)
#
# MetaboAnalystR downloads its compound and metabolite-set libraries from
# metaboanalyst.ca at run time.
# Run from repo root: Rscript src/metaboanalyst/run-untargeted-msea.R
# =============================================================================
suppressMessages({ library(dplyr); library(readr); library(stringr); library(MetaboAnalystR) })
source("src/metaboanalyst/R/run-ora.R")
source("src/metaboanalyst/R/harvest.R")

# nominal (P < 0.05) foreground, annotated with lab compound IDs and local
# mummichog candidate names
ANNOT     <- "results/milk_nominal_putative_annotation.csv"
REF_METAB <- "src/metaboanalyst/reference/reference_metabolome_1593_matched_names.txt"
OUT       <- "results/metaboanalyst/untargeted_msea"

.study_label <- function(s) dplyr::recode(s, "Misame"="MISAME-III", "Vital"="Mumta-LW",
                                          "Elicit"="ELICIT", .default = s)

# Clean an annotation string into a MetaboAnalyst-matchable compound name. Prefer
# the lab-curated compound-ID name; else the first local-mummichog putative
# candidate; else the best_annotation. Drop isotopologue / "related to" rows, strip
# the adduct suffix (_[M+H]+ ...) and any trailing "(Vitamin ...)" gloss.
.clean_name <- function(curated, putative, best) {
  nm <- ifelse(!is.na(curated) & curated != "", curated,
        ifelse(!is.na(putative) & putative != "", sub("\\s*\\|.*$", "", putative), best))
  nm <- sub("_\\[M.*$", "", nm)                 # strip adduct
  nm <- sub("\\s*\\(Vitamin[^)]*\\)", "", nm)   # strip "(Vitamin B5)" gloss
  nm <- trimws(nm)
  bad <- is.na(nm) | nm == "" | nm == "- none -" |
         str_starts(nm, "13C of") | str_starts(nm, "Related to")
  nm[bad] <- NA_character_
  nm
}

build_untargeted_query <- function() {
  d <- read.csv(ANNOT, stringsAsFactors = FALSE, check.names = FALSE)
  d$cmpd <- .clean_name(d$curated_name, d$putative_candidates, d$best_annotation)
  d %>% filter(!is.na(cmpd)) %>%
    transmute(study, studytime, visit, contrast, direction, cmpd)
}

run_untargeted_msea <- function(ref = "default", out_suffix = "", write = TRUE) {
  # ref: "default" = MetaboAnalyst whole-library background (SetMetabolomeFilter FALSE)
  #      "615"     = the 615-name reference metabolome (the published run, see below)
  #      character vector = a custom reference metabolome (sensitivity variant)
  # The background changes the p-values a lot for small query lists, so the
  # alternatives are kept for sensitivity checks.
  if (identical(ref, "default")) {
    ref_arg <- NULL
    message("background: MetaboAnalyst whole-library default (no metabolome filter)")
  } else if (identical(ref, "615")) {
    ref_arg <- trimws(readLines(REF_METAB, warn = FALSE)); ref_arg <- ref_arg[ref_arg != ""]
    message("background: 615-name reference metabolome")
  } else {
    ref_arg <- ref
    message("background: supplied reference (", length(ref), " names, suffix '", out_suffix, "')")
  }
  q <- build_untargeted_query()
  cells <- q %>% distinct(study, studytime, visit, contrast, direction)

  rows <- list()
  for (i in seq_len(nrow(cells))) {
    ce <- cells[i, ]
    query <- q %>% semi_join(ce, by = names(ce)) %>% pull(cmpd) %>% unique()
    if (length(query) < 1) next
    mSet <- tryCatch(run_ora(query, reference_names = ref_arg), error = function(e) {
      message("  skipped ", paste(unlist(ce), collapse = " / "), ": ", conditionMessage(e))
      NULL
    })
    if (is.null(mSet)) next
    res <- harvest_results(mSet)
    if (!nrow(res)) next
    res <- res %>% mutate(
      enrichment_ratio = ifelse(ce$direction == "down", -1, 1) * (hits / expected),
      study     = .study_label(ce$study),
      timepoint = ce$visit,
      contrast  = ce$contrast,
      direction = ce$direction)
    rows[[length(rows) + 1]] <- res
  }
  tab <- bind_rows(rows) %>%
    transmute(study, timepoint, contrast, direction, pathway,
              total, expected, hits, raw_p, fdr_native = fdr, enrichment_ratio) %>%
    # Pathways with < 3 members give very large enrichment ratios from one or two
    # hits. This is a reporting filter only: fdr_native is not recomputed on the
    # surviving rows (see apply_pathway_size_floor() in R/harvest.R).
    apply_pathway_size_floor(min_size = 3) %>%
    arrange(raw_p)

  if (isTRUE(write)) {
    dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
    write_csv(tab, file.path(OUT, paste0("untargeted_msea_combined", out_suffix, ".csv")))
    write_csv(filter(tab, fdr_native < 0.05),
              file.path(OUT, paste0("untargeted_msea_combined", out_suffix, "_fdr_sig.csv")))
  }
  tab
}

if (sys.nframe() == 0) {
  # The published run uses the 615-name reference metabolome as the background,
  # as described in the Methods for Fig 6A.
  tab <- run_untargeted_msea(ref = "615", write = TRUE)
  message("done untargeted MSEA: ", nrow(tab), " pathway rows; ",
          sum(tab$fdr_native < 0.05, na.rm = TRUE), " FDR-significant")
}
