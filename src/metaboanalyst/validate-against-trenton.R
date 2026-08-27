# validate-against-trenton.R — compare our scripted engine to Trenton's
# downloaded MetaboAnalyst results. Two HARD numeric gates now exist, one per
# module, both on the Elicit·Upregulated·1-month cell (the cell for which
# Trenton committed downloaded ground-truth files):
#   * validate_golden_pathway() — Pathway Analysis (pathway_results.csv), no ref
#   * validate_golden_ora()     — Enrichment/ORA (msea_ora_result.csv), with the
#                                 committed reference metabolome as background
suppressMessages({ library(dplyr); library(tibble); library(MetaboAnalystR) })
source("src/metaboanalyst/R/run-pathway.R")
source("src/metaboanalyst/R/run-ora.R")
source("src/metaboanalyst/R/harvest.R")

GOLDEN_CSV <- file.path(
  "trenton scripts/3. Results/Primary Outcomes",
  "Elicit · Upregulated · 1 month (Combined)",
  "Pathway Analysis/Download/pathway_results.csv")

# Trenton's downloaded ORA result for the SAME cell (his "Enrichment Analysis"
# upload). Vendored from the 2026 Science export alongside the pathway ground
# truth so the ORA engine has a committed numeric gate too.
GOLDEN_ORA_CSV <- file.path(
  "trenton scripts/3. Results/Primary Outcomes",
  "Elicit · Upregulated · 1 month (Combined)",
  "Enrichment Analysis (MSEA ORA)/Download/msea_ora_result.csv")

# Trenton's curated, MetaboAnalyst-ID-matched reference metabolome (the ORA
# background). Same file config-primary.R feeds the primary ORA run.
REFERENCE_PATH <- file.path(
  "trenton scripts/3. Results/Primary Outcomes",
  "Primary Outcomes Reference Metabolome.txt")

# The 7 compounds from the golden cell's downloaded datalist / Rhistory.
GOLDEN_CMPDS <- c("Nicotinic acid mononucleotide", "Niacinamide", "Nudifloramide",
                  "nicotinamide", "nicotinamide riboside",
                  "nicotinamide adenine dinucleotide", "tryptophan")

# Returns a tibble of pathways shared between our run and Trenton's downloaded
# CSV, with both raw p-values side by side.
validate_golden_pathway <- function() {
  # Our engine's pathway p-values for the golden cell's 7 compounds.
  ours <- harvest_results(run_pathway(GOLDEN_CMPDS, pathlib = "smpdb")) %>%
    transmute(pathway, raw_p_ours = raw_p)

  # Trenton's downloaded ground-truth p-values (pathways are the CSV row names).
  golden <- read.csv(GOLDEN_CSV, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  trenton <- tibble(pathway = rownames(golden), raw_p_trenton = as.numeric(golden[["Raw p"]]))

  # Keep only pathways present in both, with a relative difference for each.
  inner_join(ours, trenton, by = "pathway") %>%
    mutate(rel_diff = abs(raw_p_ours - raw_p_trenton) / pmax(raw_p_trenton, 1e-12))
}

# Returns a tibble of pathways shared between our ORA run and Trenton's
# downloaded msea_ora_result.csv, with both raw p-values side by side. Unlike the
# pathway gate, ORA restricts the background to the reference metabolome, so we
# pass the committed reference names to run_ora().
validate_golden_ora <- function() {
  reference_names <- trimws(readLines(REFERENCE_PATH, warn = FALSE))
  reference_names <- reference_names[reference_names != ""]

  # Our engine's ORA p-values for the golden cell's 7 compounds against the
  # reference-metabolome background (mset_lib "smpdb_pathway" = SMPDB 99 sets).
  ours <- harvest_results(run_ora(GOLDEN_CMPDS, reference_names = reference_names)) %>%
    transmute(pathway, raw_p_ours = raw_p)

  # Trenton's downloaded ground-truth p-values (pathways are the CSV row names).
  golden <- read.csv(GOLDEN_ORA_CSV, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  trenton <- tibble(pathway = rownames(golden), raw_p_trenton = as.numeric(golden[["Raw p"]]))

  # Keep only pathways present in both, with a relative difference for each.
  inner_join(ours, trenton, by = "pathway") %>%
    mutate(rel_diff = abs(raw_p_ours - raw_p_trenton) / pmax(raw_p_trenton, 1e-12))
}
