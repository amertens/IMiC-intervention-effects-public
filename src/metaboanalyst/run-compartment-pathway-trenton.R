# =============================================================================
# run-compartment-pathway-trenton.R
#
# The single Trenton-rule (name-based grouping) source for Table S11, Fig 6D's
# row-level data, and the pathway labels Fig 6D and Fig S10 both put on their
# y-axis. Before this, Fig 6D and Table S11/Fig S10 read DIFFERENT grouping
# rules (this script's Table S11 output vs. run-compartment-pathway.R's
# union-find `linked_upregulated_plot_data.csv`), which changed both the
# feature count and which BEP-nutrient claims the pathway table supported --
# see Manuscript/PROVENANCE_AUDIT.md §6.4. Routing all three exhibits through
# this one script's outputs closes that gap.
#
# WHY THIS EXISTS ALONGSIDE run-compartment-pathway.R
#   run-compartment-pathway.R builds a union-find `tracking_group` and names
#   each component after its lowest-P member. That is a GENERALISATION of
#   Trenton's method, not his method: it merges isobaric features he kept
#   separate, which silently drops D-lactate (m/z 135.03, absorbed by
#   threonic acid) and N1-methyl-2-pyridone-5-carboxamide (m/z 153.066,
#   absorbed by N(pi)-methyl-L-histidine) -- and with them the pyruvate,
#   glycolysis and nicotinate/nicotinamide rows the manuscript reports (the
#   niacin nutrient link). run-compartment-pathway.R and its outputs
#   (pathway_results_local.csv, linked_upregulated_plot_data.csv,
#   metabolite_pathways_local.csv) are kept ONLY as an alternate sensitivity
#   view (figS10-compartment-tracking.R's GROUPING = "union-find" switch) --
#   they no longer feed Table S11, Fig 6D, or Fig S10's primary rendering.
#
#   Trenton's rule (Compartment Tracking.Rmd 1743-1749, 2556-2562) groups by
#   the feature's OWN putative name, ported as ct_trenton_track() in
#   src/trenton-ports/compartment-tracking.R.
#
# BUILD ORDER (two-pass, see figS10-compartment-tracking.R's own header):
#   1. figure-scripts/manuscript_figures/figS10-compartment-tracking.R (1st pass)
#        -> results/compartment_tracking/trenton_pathway_compound_list.csv
#   2. this script
#        -> results/compartment_tracking/trenton_linked_crosscompartment.csv
#        -> results/compartment_tracking/metabolite_pathways_trenton.csv
#        -> results/tables/table_s11_compartment_pathway.{csv,md}
#   3. figS10-compartment-tracking.R again (2nd pass, picks up pathway labels)
#   4. figure-scripts/manuscript_figures/fig6D-crosscompartment.R
#
# Run from repo root: Rscript src/metaboanalyst/run-compartment-pathway-trenton.R
# =============================================================================
suppressMessages({ library(dplyr); library(readr); library(MetaboAnalystR) })
source("src/metaboanalyst/R/run-pathway.R")
source("src/trenton-ports/compartment-tracking.R")   # ct_trenton_track()

SUPP  <- "results/supplement_status_fdr_features.csv"
QUERY <- "results/compartment_tracking/trenton_pathway_compound_list.csv"
OUT   <- "results/compartment_tracking"
TBL   <- "results/tables"
if (!file.exists(QUERY))
  stop("run figure-scripts/manuscript_figures/figS10-compartment-tracking.R first; missing ", QUERY)

# non-informative putative name (a bare compound class, not an identity) that
# must fall back to m/z rather than be treated as a compound -- see
# run-compartment-pathway.R's identical GENERIC_NAMES / JUNK_NAMES handling.
GENERIC_NAMES <- c("Acid")

# his compartment x timepoint ordering + readable labels (Compartment
# Tracking.Rmd lines ~1893-1932); any (compartment,timepoint) not listed
# keeps its raw "comp | tp".
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

isTRUE_vec <- function(x) { if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true","1","yes") }

# ---- Fig 6D's row-level data: his rule, >= 2 distinct COMPARTMENTS ----------
# (his literal rule at 2556-2562 is >= 2 CELLS, which figS10 keeps; Fig 6D is
# specifically about cross-COMPARTMENT transfer -- the "N compartment-level
# estimates" its legend states -- so it uses the n_comp >= 2 threshold that
# ct_trenton_track() also computes.)
d0 <- read_csv(SUPP, show_col_types = FALSE) %>%
  filter(!is.na(mz), !is.na(effect_size)) %>%
  mutate(
    putative_name = { p <- trimws(putative_name); p[p %in% GENERIC_NAMES] <- ""; dplyr::na_if(p, "") },
    ct_key   = paste(compartment, timepoint, sep = " · "),
    ct_label = dplyr::coalesce(unname(CT_LEVELS[ct_key]), ct_key),
    supplement_label = case_when(
      isTRUE_vec(supplement_abundant)          ~ "Supplement-abundant",
      isTRUE_vec(supplement_reliably_detected) ~ "Detected in supplement",
      TRUE                                     ~ "Not reliably detected in supplement")) %>%
  ct_trenton_track()

linked <- d0 %>% filter(n_comp >= 2) %>%
  group_by(direction, tracking_name) %>%
  arrange(raw_p, .by_group = TRUE) %>%
  distinct(direction, tracking_name, ct_key, .keep_all = TRUE) %>%
  ungroup()

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write_csv(linked, file.path(OUT, "trenton_linked_crosscompartment.csv"))
message("cross-compartment (>=2 compartments) rows: ", nrow(linked),
        " | distinct tracking names (any direction): ", n_distinct(linked$tracking_name),
        " | upregulated tracking names: ",
        n_distinct(linked$tracking_name[linked$direction == "Upregulated"]))

# ---- KEGG pathway ORA over his literal >= 2-CELL query list -----------------
# (a superset of the >= 2-compartment set above: any compound reaching >= 2
# compartments trivially reaches >= 2 cells, so this run's per-compound hit
# membership covers every compound Fig 6D needs labels for.)
cmpds <- read_csv(QUERY, show_col_types = FALSE)$compound
message("Table S11 query compounds (>=2 cells, his literal rule): ", length(cmpds))

h <- harvest_pathway(run_pathway(cmpds, pathlib = "kegg"))
dir.create(TBL, recursive = TRUE, showWarnings = FALSE)
write_csv(h$metabolite_pathways, file.path(OUT, "metabolite_pathways_trenton.csv"))

# The driving-metabolite column IS the over-representation hit membership of
# this same run, so the statistics and the driver list cannot disagree (the
# printed Table S11 took that column from a hand-coded map that did not match
# its own run -- Manuscript/PROVENANCE_AUDIT.md §1.1a).
drivers <- h$metabolite_pathways %>% group_by(kegg_id) %>%
  summarise(driving = paste(sort(unique(tracking_name)), collapse = "; "), .groups = "drop")
s11 <- h$pathways %>% left_join(drivers, by = "kegg_id") %>% arrange(raw_p) %>%
  transmute(pathway, kegg_id, total, expected = round(expected, 3), hits,
            raw_p = signif(raw_p, 3), fdr = signif(fdr, 3), impact = round(impact, 3),
            driving_putative_metabolites = driving)
write_csv(s11, file.path(TBL, "table_s11_compartment_pathway.csv"))

fmt_p <- function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
writeLines(c(
  "| Pathway | Total | Expected | Hits | *P* | FDR | Impact | Driving putative metabolites |",
  "|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %d | %.3f | %d | %s | %s | %s | %s |",
          s11$pathway, s11$total, s11$expected, s11$hits, fmt_p(s11$raw_p), fmt_p(s11$fdr),
          formatC(s11$impact, format = "g", digits = 3),
          ifelse(is.na(s11$driving_putative_metabolites), "—", s11$driving_putative_metabolites))),
  file.path(TBL, "table_s11_compartment_pathway.md"))

cat("query q (mapped) =", round(h$pathways$expected[1] / h$pathways$total[1] * 1592, 2),
    "| pathways:", nrow(s11), "\n-> ", file.path(TBL, "table_s11_compartment_pathway.csv"), "\n")
