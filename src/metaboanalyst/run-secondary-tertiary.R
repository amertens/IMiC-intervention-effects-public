# run-secondary-tertiary.R, populate the secondary + tertiary outcome-group
# directories (combined and stratified arms).
#
# MODULE PER GROUP (must match how each group is REPORTED):
#   secondary -> Pathway module (pathora, filter OFF).
#   tertiary  -> ORA module (msetora) with the WHOLE-METABOLOME reference (615
#                library names), IDENTICAL to the Fig 5B tertiary MSEA
#                (run-tertiary-msea.R). This was previously (incorrectly) run
#                through the pathway module with a filter-OFF whole-library
#                background, so the tertiary "ORA by direction" supplement (Fig S8)
#                and its supplementary table disagreed with Fig 5B on both the test
#                (pathway impact vs ORA) and the universe. Running tertiary as ORA
#                here makes results/metaboanalyst/tertiary_{combined,stratified}/
#                and everything derived from them (the supplementary table, the
#                direction-split figure, the replication matrix) consistent with 5B.
#
# Writes results/metaboanalyst/<group>_<arm>/ (per-cell + combined + skip log).
suppressMessages({ library(dplyr) })
source("src/metaboanalyst/R/run-outcome-group.R")

COMBINED_RDS   <- "results/combined_intervention_effects_results_combined_arms.RDS"
STRATIFIED_RDS <- "results/combined_intervention_effects_results_stratified_arms.RDS"
# Same reference metabolome (1,593 HMDB -> 615 library names) used by Fig 5B.
REF_METAB      <- "src/metaboanalyst/reference/reference_metabolome_1593_matched_names.txt"

.reference_metabolome <- function() {
  if (!file.exists(REF_METAB)) {
    stop("Reference metabolome not found at '", REF_METAB, "'. ",
         "Build it with src/metaboanalyst/build-reference-metabolome.R.",
         call. = FALSE)
  }
  ref <- trimws(readLines(REF_METAB, warn = FALSE))
  ref[ref != ""]
}

run_secondary_tertiary <- function(write = TRUE, max_cells = Inf) {
  combined   <- readRDS(COMBINED_RDS)
  stratified <- readRDS(STRATIFIED_RDS)
  ref_metabolome <- .reference_metabolome()
  message("tertiary ORA background: reference metabolome, ", length(ref_metabolome), " names")

  runs <- list()

  # --- secondary: Pathway module (unchanged) ---------------------------------
  message("== secondary / combined ==")
  runs[["secondary_combined"]] <-
    run_outcome_group(combined, "secondary", module = "pathway", arm_set = "combined",
                      write = write, max_cells = max_cells)
  message("== secondary / stratified ==")
  runs[["secondary_stratified"]] <-
    run_outcome_group(stratified, "secondary", module = "pathway", arm_set = "stratified",
                      write = write, max_cells = max_cells)

  # --- tertiary: ORA module + reference metabolome (matches Fig 5B) -----
  message("== tertiary / combined (ORA + 615 background) ==")
  runs[["tertiary_combined"]] <-
    run_outcome_group(combined, "tertiary", module = "ora", arm_set = "combined",
                      reference_names = ref_metabolome, write = write, max_cells = max_cells)
  message("== tertiary / stratified (ORA + 615 background) ==")
  runs[["tertiary_stratified"]] <-
    run_outcome_group(stratified, "tertiary", module = "ora", arm_set = "stratified",
                      reference_names = ref_metabolome, write = write, max_cells = max_cells)

  # brief console summary (guard groups that produced no results: run_outcome_group
  # returns a 0-column tibble() when every cell is skipped, so only subset when non-empty)
  for (run_name in names(runs)) {
    run <- runs[[run_name]]
    n_cells_with_results <- if (nrow(run$results) > 0)
      nrow(dplyr::distinct(run$results[c("study", "contrast", "direction")])) else 0L
    message(sprintf("  %-24s cells-with-results=%d  skipped=%d",
                    run_name, n_cells_with_results, nrow(run$skipped)))
  }
  runs
}

if (sys.nframe() == 0) invisible(run_secondary_tertiary(write = TRUE))
