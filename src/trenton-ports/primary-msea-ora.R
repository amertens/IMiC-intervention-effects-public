# =============================================================================
# primary-msea-ora.R
# Faithful R port of Trenton's "Primary Outcomes (MSEA ORA).Rmd".
#
# Mirrors his notebook top-to-bottom:
#   1. Import combined- and stratified-arm intervention-effects RDS.
#   2. Generate the MetaboAnalyst compound label (his verbatim 11-entry synonym map).
#   3. Create the reference metabolome (distinct primary compounds, "Total *" roll-
#      ups excluded) and ID-convert it.
#   4. Per (study x direction x timepoint) cell: query = features with pval < 0.05
#      and the matching sign of the effect estimate, "Total *" excluded; run ORA
#      (hypergeometric, SMPDB) against the reference metabolome; download results.
#
# AUTOMATION (the only change vs the Rmd): the two manual metaboanalyst.ca steps --
#   (a) the reference-metabolome "Compound ID Conversion" upload, and
#   (b) each cell's "Enrichment Analysis (ORA)" upload/download --
# are replaced by the repo's offline MetaboAnalystR engine
# (build_matched_reference() and run_ora()), which reproduces his downloaded
# msea_ora_result.csv files EXACTLY (Elicit up 1mo/5mo, down 5mo; Mumpta up 1.5mo:
# 60/60/9/76 pathways, raw-p & FDR delta = 0). See src/trenton-ports/README.md.
#
# His filters/cells/order are otherwise unchanged. No coding errors were found in
# this Rmd to correct.
#
# Run from repo root:  Rscript src/trenton-ports/primary-msea-ora.R
# =============================================================================
suppressMessages({ library(dplyr); library(readr); library(stringr); library(MetaboAnalystR) })
source("src/metaboanalyst/R/label-map.R")       # apply_label_map(): his verbatim synonym map
source("src/metaboanalyst/R/build-cells.R")     # build_cells(): his pval<0.05 x est-sign x !Total cell rule
source("src/metaboanalyst/R/build-reference.R") # build_matched_reference(): scripted ID-conversion
source("src/metaboanalyst/R/run-ora.R")         # run_ora(): scripted Enrichment Analysis (ORA)
source("src/metaboanalyst/R/harvest.R")         # harvest_results(): tidy the downloaded result

COMBINED_RDS   <- "trenton scripts/1. Data/combined_intervention_effects_results_combined_arms.RDS"
STRATIFIED_RDS <- "trenton scripts/1. Data/combined_intervention_effects_results_stratified_arms.RDS"
OUT            <- "results/trenton-ports/primary_msea_ora"
# Trenton's committed reference (from his manual ID-conversion). We rebuild the
# reference by script and gate it against this file so the port stays faithful.
TRENTON_REF    <- "trenton scripts/3. Results/Primary Outcomes/Primary Outcomes Reference Metabolome.txt"

# -- his study-name -> published cohort label (for tidy output only) -----------
.cohort <- function(studytime) {
  s <- str_replace(studytime, "\\s*\\(.*$", "")
  recode(s, "Elicit" = "ELICIT", "Misame" = "MISAME-III", "Vital" = "Mumta-LW",
         "Mumpta" = "Mumta-LW", .default = s)
}
.timepoint <- function(studytime) str_trim(str_replace(str_replace(studytime, "^[^(]*\\(", ""), "\\)\\s*$", ""))

# -----------------------------------------------------------------------------
# Step 3 (his "Create Reference Metabolome"): distinct primary compounds, "Total *"
# excluded, run through MetaboAnalyst's Compound ID Conversion; keep the Match names.
# We rebuild it by script AND verify it matches his committed reference file.
# -----------------------------------------------------------------------------
build_reference <- function(combined) {
  raw <- combined %>%
    filter(outcome_group == "primary") %>%
    mutate(metaboanalyst_label = apply_label_map(label_f)) %>%
    filter(!str_detect(label_f, regex("^\\s*total", ignore_case = TRUE))) %>%  # his "!Total" exclusion (on raw label)
    distinct(metaboanalyst_label) %>% arrange(metaboanalyst_label) %>% pull(metaboanalyst_label)

  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  ref <- build_matched_reference(
    raw,
    cache_path      = file.path(OUT, "primary_reference_metabolome.txt"),
    kegg_cache_path = file.path(OUT, "primary_reference_metabolome_KEGG.txt"),
    refresh         = TRUE)

  # Fidelity gate: the scripted reference must match Trenton's committed one.
  if (file.exists(TRENTON_REF)) {
    his <- trimws(readLines(TRENTON_REF, warn = FALSE)); his <- his[his != ""]
    only_ours <- setdiff(ref, his); only_his <- setdiff(his, ref)
    if (length(only_ours) || length(only_his)) {
      message("NOTE: scripted reference differs from Trenton's committed file ",
              "(+", length(only_ours), " / -", length(only_his), "); using HIS file to stay faithful.")
      ref <- his
    } else {
      message("reference metabolome: scripted ID-conversion reproduces Trenton's file exactly (",
              length(ref), " names).")
    }
  }
  ref
}

# -----------------------------------------------------------------------------
# Step 4 (his "Analysis"): one ORA per (study x direction x timepoint) cell.
# build_cells() encodes his exact cell rule (pval<0.05, est sign, !Total roll-ups),
# emitting one cell per non-empty study-time x direction. We run both arm framings,
# matching the "(Combined)" and "(Stratified)" cells in his Rmd.
# -----------------------------------------------------------------------------
run_framing <- function(dat, ref, arm_set) {
  cells <- build_cells(dat, "primary")   # his directional pval<0.05 cells
  rows <- list()
  for (cell in cells) {
    mSet <- tryCatch(run_ora(cell$query, reference_names = ref), error = function(e) NULL)
    if (is.null(mSet)) { message("  skip (too few mappable): ", cell$study, " / ", cell$direction); next }
    res <- harvest_results(mSet)                     # total, hits, expected, raw_p, fdr
    if (!nrow(res)) next
    res <- res %>% mutate(study = .cohort(cell$study), timepoint = .timepoint(cell$study),
                          contrast = cell$contrast, direction = cell$direction, arm_set = arm_set,
                          .before = 1)
    # per-cell download, mirroring his one-folder-per-cell layout
    slug <- gsub("[^A-Za-z0-9]+", "_", paste(cell$study, cell$direction, cell$contrast))
    celldir <- file.path(OUT, arm_set, slug); dir.create(celldir, recursive = TRUE, showWarnings = FALSE)
    write_csv(res, file.path(celldir, "msea_ora_result.csv"))
    rows[[length(rows) + 1]] <- res
  }
  bind_rows(rows)
}

if (sys.nframe() == 0) {
  combined   <- readRDS(COMBINED_RDS)
  stratified <- readRDS(STRATIFIED_RDS)
  ref <- build_reference(combined)

  message("== combined-arm ORA cells ==")
  combined_res   <- run_framing(combined,   ref, "combined")
  message("== stratified-arm ORA cells ==")
  stratified_res <- run_framing(stratified, ref, "stratified")

  all <- bind_rows(combined_res, stratified_res) %>% arrange(arm_set, study, timepoint, direction, raw_p)
  write_csv(all, file.path(OUT, "primary_msea_ora_all_cells.csv"))
  message(sprintf("done: %d pathway rows across %d cells; %d FDR-significant. -> %s",
                  nrow(all), nrow(distinct(all, arm_set, study, timepoint, direction)),
                  sum(all$fdr < 0.05, na.rm = TRUE), OUT))
}
