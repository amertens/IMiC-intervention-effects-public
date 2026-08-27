# run-compartment-pathway.R -- SUPERSEDED / SENSITIVITY-VIEW ONLY.
#
# This union-find grouping was a GENERALISATION of Trenton's published method,
# not his method (Manuscript/PROVENANCE_AUDIT.md §6). It no longer feeds Fig 6D,
# Table S11, or Fig S10's primary rendering -- all three now read the name-based
# rule built by src/metaboanalyst/run-compartment-pathway-trenton.R (Fig 6D,
# Table S11) and src/trenton-ports/compartment-tracking.R's ct_trenton_track()
# (Fig S10). This script and its outputs (pathway_results_local.csv,
# linked_upregulated_plot_data.csv, metabolite_pathways_local.csv) are kept
# ONLY so figS10-compartment-tracking.R's GROUPING = "union-find" switch has
# something to render for comparison. Do not point any manuscript exhibit at
# this script's outputs.
#
# Historical header, for that comparison view -- web-independent pathway step
# for Figure 6D (cross-compartment tracking of BEP-associated UPREGULATED
# metabolite features).
#
# Replaces the TWO web/manual inputs in Trenton's "Compartment Tracking.Rmd":
#   (1) the downloaded MetaboAnalyst "pathway_results.csv"      (## Figure With Pathway Results)
#   (2) the hardcoded metabolite->pathway `tribble`             (metaboanalyst_metabolite_pathways)
# ...with a LOCAL MetaboAnalystR KEGG pathway analysis (src/metaboanalyst/R/run-pathway.R)
# on exactly his compound list, and the ORA hit-membership as the metabolite->pathway map.
#
# Faithful to his linkage (Compartment Tracking.Rmd, `compartment_pairs` +
# `linked_compartment_tracking_plot_data`):
#   - features = the FDR-significant cross-compartment set with supplement status
#     (results/supplement_status_fdr_features.csv, from src/2 analysis/54-supplement-detection-status.R),
#   - features are grouped into "the same metabolite" by his THREE-TIER same-direction
#     linkage (exact feature id / identical normalized putative name / same ion mode &
#     mass within 25 ppm), ported in src/trenton-ports/compartment-tracking.R
#     (ct_tracking_components); the group's display tracking_name = a real putative name
#     if any group member has one, else "m/z <round4>",
#   - keep groups significant in >= 2 distinct COMPARTMENTS (his concordant cross-
#     compartment transfer; NOT merely >= 2 compartment x timepoint cells, which the
#     earlier name/m/z-string recurrence rule allowed to be met within one compartment),
#   - the pathway compound list = distinct UPREGULATED group names that are named
#     (not "m/z ...") -- his `metaboanalyst_pathway_compond_list`.
#
# Outputs (all in-repo, no web):
#   results/compartment_tracking/linked_upregulated_plot_data.csv   (per compound x compartment-timepoint, for the figure)
#   results/compartment_tracking/pathway_results_local.csv          (replaces web pathway_results.csv)
#   results/compartment_tracking/metabolite_pathways_local.csv      (replaces the tribble)
#
# Run from repo root:  Rscript src/metaboanalyst/run-compartment-pathway.R
suppressMessages({ library(dplyr); library(stringr); library(readr); library(tidyr); library(MetaboAnalystR) })
source("src/metaboanalyst/R/run-pathway.R")
# Trenton's three-tier same-direction cross-compartment linkage (exact feature-id /
# identical normalized putative name / same ion mode & mass within 25 ppm), ported in
# src/trenton-ports/compartment-tracking.R. Sourced for ct_tracking_components(). NOTE it
# defines its own OUT/SLIM/COMP_LEVELS, so we (re)define our paths AFTER this source.
source("src/trenton-ports/compartment-tracking.R")

SUPP <- "results/supplement_status_fdr_features.csv"
OUT  <- "results/compartment_tracking"

# Trenton's compartment x timepoint ordering + readable labels (Compartment Tracking.Rmd
# lines ~1893-1932). Any (compartment,timepoint) not listed keeps its raw "comp | tp".
CT_LEVELS <- c(
  "Maternal plasma · incl" = "Maternal plasma · Enrollment",
  "Maternal plasma · tri3" = "Maternal plasma · Third trimester",
  "Maternal plasma · pn12" = "Maternal plasma · 1-2 months",
  "Maternal VAMS · tri3"   = "Maternal VAMS · Third trimester",
  "Maternal VAMS · pn56"   = "Maternal VAMS · 5-6 months",
  "Milk · 1421d"           = "Milk · 14-21 days",
  "Milk · pn12"            = "Milk · 1-2 months",
  "Milk · pn34"            = "Milk · 3-4 months",
  "Infant VAMS · acco"     = "Infant VAMS · Birth",
  "Infant VAMS · pn12"     = "Infant VAMS · 1-2 months",
  "Infant VAMS · pn34"     = "Infant VAMS · 3-4 months",
  "Infant VAMS · pn56"     = "Infant VAMS · 5-6 months")

# robust logical coercion (the flags read back as logical, but guard NA/character)
isTRUE_vec <- function(x) { if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true","1","yes") }

# Non-informative putative names (a bare compound class, not an identity) that must
# fall back to m/z rather than be treated as a compound. "Acid" appears as a stray
# putative_name in the FDR set; without this it would label a row and be sent to KEGG.
GENERIC_NAMES <- c("Acid")

build_linked <- function() {
  d0 <- read_csv(SUPP, show_col_types = FALSE) %>%
    filter(!is.na(mz), !is.na(effect_size)) %>%
    mutate(
      pn0 = { p <- trimws(putative_name); p[p %in% GENERIC_NAMES] <- ""; dplyr::na_if(p, "") },
      ct_key   = paste(compartment, timepoint, sep = " · "),
      ct_label = dplyr::coalesce(unname(CT_LEVELS[ct_key]), ct_key),
      supplement_label = case_when(
        isTRUE_vec(supplement_abundant)          ~ "Supplement-abundant",
        isTRUE_vec(supplement_reliably_detected) ~ "Detected in supplement",
        TRUE                                     ~ "Not reliably detected in supplement"))

  # FAITHFUL LINKAGE. SUPP is already the FDR-significant feature set, so pass every row
  # through Trenton's three-tier same-direction union-find; ct_tracking_components() adds
  # `tracking_group` (a connected component = "one metabolite"). This replaces the earlier
  # putative-name / m/z-string recurrence, which had no feature-id or 25-ppm mass tier and
  # could split concordant features labelled only by (divergent) m/z.
  d0 <- ct_tracking_components(d0)

  # One display name per tracking_group: a real putative name if ANY group member has one
  # (taken from the strongest / lowest-raw_p member), else the group's m/z. Shared across
  # the group so its features collapse to ONE figure row and ONE KEGG compound.
  grp_name <- d0 %>% group_by(tracking_group) %>%
    arrange(raw_p, .by_group = TRUE) %>%
    summarise(pn  = { x <- pn0[!is.na(pn0)]; if (length(x)) x[1] else NA_character_ },
              mz0 = dplyr::first(mz), .groups = "drop") %>%
    mutate(tracking_name = dplyr::coalesce(pn, paste0("m/z ", round(mz0, 4)))) %>%
    select(tracking_group, tracking_name)

  d <- d0 %>% left_join(grp_name, by = "tracking_group")

  # Concordant cross-compartment transfer: keep groups significant in >= 2 distinct
  # COMPARTMENTS. One row per (group, compartment x timepoint cell), strongest = min raw_p.
  keep <- d %>% distinct(tracking_group, compartment) %>%
    count(tracking_group) %>% filter(n >= 2) %>% pull(tracking_group)

  d %>% filter(tracking_group %in% keep) %>%
    group_by(direction, tracking_group) %>%
    arrange(raw_p, .by_group = TRUE) %>%
    distinct(direction, tracking_group, ct_key, .keep_all = TRUE) %>%
    ungroup()
}

if (sys.nframe() == 0) {
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  linked <- build_linked()
  write_csv(linked, file.path(OUT, "linked_upregulated_plot_data.csv"))  # both directions kept; figure filters up

  cmpds <- linked %>%
    filter(direction == "Upregulated", !str_starts(tracking_name, "m/z")) %>%
    distinct(tracking_name) %>% pull(tracking_name) %>% sort()
  message("upregulated named compounds for pathway analysis: ", length(cmpds))

  mSet <- run_pathway(cmpds, pathlib = "kegg")
  h <- harvest_pathway(mSet)
  write_csv(h$pathways, file.path(OUT, "pathway_results_local.csv"))
  write_csv(h$metabolite_pathways, file.path(OUT, "metabolite_pathways_local.csv"))

  sig <- h$pathways %>% filter(raw_p < 0.05)
  message("pathways: ", nrow(h$pathways), " | raw p<0.05: ", nrow(sig),
          " | metabolite-pathway links: ", nrow(h$metabolite_pathways))

  # ---- Table S11 (printed supplement) -------------------------------------
  # The printed table pairs each pathway's ORA statistics with the metabolites
  # DRIVING it. Trenton's version took that last column from a HAND-CODED map that
  # did not match the run producing the P-values (it credited glyoxylate to
  # "Citrate; Mesaconate" although the hypergeometric P of 0.019 could only have
  # come from a run in which Mesaconate was unmapped -- see
  # Manuscript/PROVENANCE_AUDIT.md 1.1). Here the column IS the ORA hit membership,
  # so the statistics and the driver list can never disagree again.
  TBL <- "results/tables"
  dir.create(TBL, recursive = TRUE, showWarnings = FALSE)
  drivers <- h$metabolite_pathways %>% group_by(kegg_id) %>%
    summarise(driving = paste(sort(unique(tracking_name)), collapse = "; "), .groups = "drop")
  s11 <- h$pathways %>% left_join(drivers, by = "kegg_id") %>%
    arrange(raw_p) %>%
    transmute(pathway, kegg_id, total, expected = round(expected, 3), hits,
              raw_p = signif(raw_p, 3), fdr = signif(fdr, 3),
              impact = round(impact, 3), driving_putative_metabolites = driving)
  write_csv(s11, file.path(TBL, "table_s11_compartment_pathway.csv"))
  fmt_p <- function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  writeLines(c(
    "| Pathway | Total | Expected | Hits | *P* | FDR | Impact | Driving putative metabolites |",
    "|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %d | %.3f | %d | %s | %s | %s | %s |",
            s11$pathway, s11$total, s11$expected, s11$hits, fmt_p(s11$raw_p),
            fmt_p(s11$fdr), formatC(s11$impact, format = "g", digits = 3),
            ifelse(is.na(s11$driving_putative_metabolites), "—",
                   s11$driving_putative_metabolites))),
    file.path(TBL, "table_s11_compartment_pathway.md"))
  message("-> ", OUT, " ; ", TBL, "/table_s11_compartment_pathway.{csv,md}")
}
