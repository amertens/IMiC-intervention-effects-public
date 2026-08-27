# run-untargeted-msea.R -- reproduce the manuscript's untargeted milk MSEA (Fig 6A)
# fully in-repo with NO web tool (2026-08-13 web-independence restore).
#
# Pipeline (all local):
#   raw untargeted ATE (results/adjusted_*_untargeted_*_ATE.RDS)
#     -> src/2 analysis/47-annotate-milk-features.R  [lab compound-ID files +
#        LOCAL mummichog candidate names (results/mummichog_output_directional/),
#        NOT the MetaboAnalyst web matcher]
#     -> results/milk_fdr_sig_putative_annotation.csv  (query = named FDR-sig features)
#   this script: MetaboAnalystR ORA (offline) of each cell's named query against the
#     REFERENCE METABOLOME background (615 library-matched names; the same background
#     that reproduces the submitted signals, e.g. spermidine Q = 2.3e-6).
#
# Background note: the earlier version of this script read Trenton's DOWNLOADED
# MetaboAnalyst web ORA tables from data/msea/ (a measured-panel background under
# which nothing clears FDR). This restores the local, reference-background ORA.
# Run from repo root:  Rscript src/metaboanalyst/run-untargeted-msea.R
suppressMessages({ library(dplyr); library(readr); library(stringr); library(MetaboAnalystR) })
source("src/metaboanalyst/R/run-ora.R")
source("src/metaboanalyst/R/harvest.R")

# NOMINAL (P<0.05) foreground -- the ORA foreground the submitted Fig 6A method uses
# ("untargeted features with unadjusted P < 0.05"). Built locally by
# src/2 analysis/47-annotate-milk-features.R (lab compound-IDs + local mummichog).
ANNOT     <- "results/milk_nominal_putative_annotation.csv"
# Background = MetaboAnalyst WHOLE-LIBRARY default (SetMetabolomeFilter FALSE); see the
# 2026-08-14 reconciliation below. The DEFAULT is intentionally the non-inflated universe.
#
# 2026-08-14 RECONCILIATION (against Trenton's actual submitted web-ORA outputs
# msea_ora_result_*.xlsx, recovered from Andrew's Downloads):
#   - The submitted Fig 6A = per-cell SMPDB ("pathway-associated metabolite sets") ORA on
#     the nominal-P<0.05 foreground (sig==1, CONFIRMED identical to unadjusted p<0.05).
#     Same library + same threshold as this script.
#   - The excess significance of the earlier 615-name run is NOT the method: it is the
#     ANNOTATION + REFERENCE-METABOLOME layer. The 615 reference under/over-covers specific
#     pathways vs Trenton's hand-annotation, inflating a few cells (Bile Acid Biosynthesis:
#     submitted 2 hits/total 5 -> p=0.47 n.s.; the 615 run gave 5 hits/total 15 -> p=7e-5).
#     Background SWEEP (foreground fixed): whole-library default best p~2e-4 / ER~7 (matches
#     the submitted panel's modest magnitude); 615 -> p~1e-6 / ER~47; QER-1270 -> p~3e-9.
#   - So the DEFAULT here is the whole-library universe: reproducible, un-inflated, and the
#     closest reproducible match to the submitted panel. These are EXPLORATORY second-platform
#     signals -- read the pathway THEMES/directions, not the exact p. The submission's exact
#     numbers depend on hand-annotation keys not in the repo and are not recoverable.
#   Pass ref = "615" (legacy small-N reference) or a custom character vector for a
#   sensitivity variant; ref = "default" (the default) = whole-library.
REF_METAB <- "src/metaboanalyst/reference/reference_metabolome_1593_matched_names.txt"  # legacy 615-name (ref="615")
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
  # ref: "default" = MetaboAnalyst WHOLE-LIBRARY background (SetMetabolomeFilter FALSE) --
  #        the un-inflated, submission-faithful universe (2026-08-14 reconciliation).
  #      "615"     = legacy 615-name reference metabolome (small-N; inflates significance).
  #      character vector = a custom reference metabolome (sensitivity variant).
  if (identical(ref, "default")) {
    ref_arg <- NULL
    message("background: MetaboAnalyst whole-library default (no metabolome filter)")
  } else if (identical(ref, "615")) {
    ref_arg <- trimws(readLines(REF_METAB, warn = FALSE)); ref_arg <- ref_arg[ref_arg != ""]
    message("background: legacy 615-name reference metabolome")
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
    mSet <- tryCatch(run_ora(query, reference_names = ref_arg), error = function(e) NULL)
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
    # 2026-08-26 (Option C, author-approved): >=3-member pathway size floor, matching
    # the convention already used for the tertiary (5B) and supplementary ORA tables
    # (Table S4/S6) elsewhere in this repo. This is a REPORTING filter only -- fdr_native
    # is MetaboAnalyst's per-cell BH FDR computed over the FULL reference-restricted
    # library before this filter, and is NOT recomputed on the surviving rows (see
    # apply_pathway_size_floor()'s own header: recomputing on the filtered set would
    # shrink the testing family and inflate significance). Previously min_size = 1 (no
    # floor), which let single-digit "total" pathways (e.g. 2/2, 3/3 hit fractions)
    # dominate the labelled/tabulated results -- a small-N ORA artifact, not signal.
    apply_pathway_size_floor(min_size = 3) %>%
    arrange(raw_p)

  if (isTRUE(write)) {
    dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
    write_csv(tab, file.path(OUT, paste0("untargeted_msea_combined", out_suffix, ".csv")))
    write_csv(filter(tab, fdr_native < 0.05),
              file.path(OUT, paste0("untargeted_msea_combined", out_suffix, "_perContrastFDRsig_tableS4.csv")))
  }
  tab
}

if (sys.nframe() == 0) {
  # 2026-08-26 (Option C, author-approved): ref = "615" (the constructed reference
  # metabolome, 1,593 HMDB-matched measured features -> 615 library-mappable names).
  # This is what the manuscript Methods already describes ("we constructed a reference
  # metabolome by integrating study-specific metabolite feature labels... untargeted-
  # metabolome enrichment (Fig. 6A) used the reference metabolome described above as
  # the enrichment background") -- the previous ref = "default" (whole-library) run
  # did NOT match that description, and additionally zeroed out MISAME-III entirely
  # (0 FDR-significant pathways vs 8/8 for ELICIT/Mumta-LW), which the 615 reference
  # does not do (see Manuscript/RESUBMISSION_PACKAGE_FIG6D_AND_TODOS_2026-08-25.md).
  tab <- run_untargeted_msea(ref = "615", write = TRUE)
  message("done untargeted MSEA: ", nrow(tab), " pathway rows; ",
          sum(tab$fdr_native < 0.05, na.rm = TRUE), " FDR-significant")
}
