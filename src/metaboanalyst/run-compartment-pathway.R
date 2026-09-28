# =============================================================================
# run-compartment-pathway.R
#
# Cross-compartment pathway analysis for MISAME-III (Table S11) and the row data
# and pathway labels for Fig 6D. Groups the FDR-significant features by their own
# putative name (ct_name_track()) and writes the features whose name reaches at
# least two compartments (Fig 6D rows). It then runs a KEGG pathway analysis
# (MetaboAnalystR, hypergeometric test) on the query list from
# build-compartment-query-list.R and writes Table S11, whose driving-metabolite
# column is taken from the same run's hit membership.
#
# Inputs : results/supplement_status_fdr_features.csv  (src/2 analysis/54-supplement-detection-status.R)
#          results/compartment_tracking/pathway_compound_list.csv  (build-compartment-query-list.R)
#          src/metaboanalyst/reference/kegg_hsa_pathway_names.csv  (via R/run-pathway.R)
# Outputs: results/compartment_tracking/linked_crosscompartment.csv  (Fig 6D rows)
#          results/compartment_tracking/metabolite_pathways.csv      (compound -> pathway hits; Fig 6D labels)
#          results/tables/table_s11_compartment_pathway.csv          (Table S11)
#          results/tables/table_s11_compartment_pathway.md           (same table, Markdown)
#
# Build order: build-compartment-query-list.R -> this script ->
#   figure-scripts/manuscript_figures/fig6D-crosscompartment.R
# MetaboAnalystR downloads its KEGG library from metaboanalyst.ca at run time.
# Run from repo root: Rscript src/metaboanalyst/run-compartment-pathway.R
# =============================================================================
suppressMessages({ library(dplyr); library(readr); library(MetaboAnalystR) })
source("src/metaboanalyst/R/run-pathway.R")
source("src/metaboanalyst/R/compartment-tracking.R")   # ct_name_track()

SUPP  <- "results/supplement_status_fdr_features.csv"
QUERY <- "results/compartment_tracking/pathway_compound_list.csv"
OUT   <- "results/compartment_tracking"
TBL   <- "results/tables"
if (!file.exists(QUERY))
  stop("run src/metaboanalyst/build-compartment-query-list.R first; missing ", QUERY)

# A bare compound class, not an identity: such features fall back to their m/z
# label rather than being grouped as one compound (the query list drops the same
# name).
GENERIC_NAMES <- c("Acid")

# compartment x timepoint order and readable labels from the original R Markdown
# analysis; any (compartment, timepoint) not listed keeps its raw "comp · tp".
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

# ---- Fig 6D row data: names reaching >= 2 distinct compartments -------------
# The query list below uses >= 2 compartment x timepoint cells; Fig 6D shows
# transfer across compartments, so it uses the n_comp >= 2 count that
# ct_name_track() also computes. Keep the lowest-P feature per name x cell.
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
  ct_name_track()

linked <- d0 %>% filter(n_comp >= 2) %>%
  group_by(direction, tracking_name) %>%
  arrange(raw_p, .by_group = TRUE) %>%
  distinct(direction, tracking_name, ct_key, .keep_all = TRUE) %>%
  ungroup()

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write_csv(linked, file.path(OUT, "linked_crosscompartment.csv"))
message("cross-compartment (>=2 compartments) rows: ", nrow(linked),
        " | distinct tracking names (any direction): ", n_distinct(linked$tracking_name),
        " | upregulated tracking names: ",
        n_distinct(linked$tracking_name[linked$direction == "Upregulated"]))

# ---- KEGG pathway analysis over the >= 2-cell query list --------------------
# Any name reaching >= 2 compartments also reaches >= 2 cells, so this run's hit
# membership covers every compound Fig 6D needs a pathway label for.
cmpds <- read_csv(QUERY, show_col_types = FALSE)$compound
message("Table S11 query compounds (>= 2 compartment x timepoint cells): ", length(cmpds))

h <- harvest_pathway(run_pathway(cmpds, pathlib = "kegg"))
dir.create(TBL, recursive = TRUE, showWarnings = FALSE)
write_csv(h$metabolite_pathways, file.path(OUT, "metabolite_pathways.csv"))

# The driving-metabolite column is the hit membership of this same run, so the
# statistics and the driver list cannot disagree.
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
