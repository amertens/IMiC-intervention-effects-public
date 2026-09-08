# run-outcome-group.R, general driver: run one outcome group through one module.
#
# Dispatches ORA (msetora, with an optional reference metabolome) or Pathway
# Analysis (pathora, filter OFF). Per Trenton's rule: primary -> ORA,
# secondary/tertiary/exploratory -> pathway. Thin cells (pathway needs >= 3
# mapped metabolites) are caught and logged, not fatal.
suppressMessages({ library(dplyr); library(MetaboAnalystR) })
source("src/metaboanalyst/R/build-cells.R")
source("src/metaboanalyst/R/run-ora.R")
source("src/metaboanalyst/R/run-pathway.R")
source("src/metaboanalyst/R/harvest.R")

# Returns list(results = <tibble>, skipped = <tibble>, out_dir = <path>).
run_outcome_group <- function(dat, outcome_group,
                              module = c("ora", "pathway"),
                              arm_set = "combined",
                              reference_names = NULL,
                              query_mode = c("directional_pval", "combined_sigfdr"),
                              write = TRUE, max_cells = Inf,
                              results_root = "results/metaboanalyst") {
  module <- match.arg(module)
  query_mode <- match.arg(query_mode)
  # "directional_pval" (default): directional pval<0.05 cells (ORA / MSEA family).
  # "combined_sigfdr": one non-directional sigFDR==1 cell per study x contrast
  #   (Trenton's pathway-impact construction; used for Fig 3B).
  cells <- if (query_mode == "combined_sigfdr")
             build_cells_combined_sigfdr(dat, outcome_group)
           else
             build_cells(dat, outcome_group)
  if (is.finite(max_cells)) cells <- head(cells, max_cells)
  out_dir <- file.path(results_root, paste0(outcome_group, "_", arm_set))

  results_list <- list()
  skipped <- list()

  for (cell in cells) {
    # Filesystem-safe folder name for this cell's per-cell outputs, e.g.
    # "Elicit_1_mo_up_Nico". Any non-alphanumeric run collapses to a single "_".
    slug <- gsub("[^A-Za-z0-9]+", "_", paste(cell$study, cell$direction, cell$contrast))
    # A cell can fail legitimately (too few mappable metabolites); we catch that
    # below and record it in `skipped` rather than aborting the whole group.
    tryCatch({
      mSet <- if (module == "ora") run_ora(cell$query, reference_names = reference_names)
              else                 run_pathway(cell$query)

      # Tag every output row with the cell's identity so results from all cells
      # can be row-bound and still be traceable back to their study/contrast.
      cell_results <- harvest_results(mSet) %>%
        mutate(study = cell$study, contrast = cell$contrast, direction = cell$direction, .before = 1)
      cell_membership <- harvest_membership(mSet) %>%
        mutate(study = cell$study, contrast = cell$contrast, direction = cell$direction, .before = 1)

      if (isTRUE(write)) {
        cell_dir <- file.path(out_dir, slug)
        dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
        write.csv(cell_results, file.path(cell_dir, "results.csv"), row.names = FALSE)
        write.csv(cell_membership, file.path(cell_dir, "membership.csv"), row.names = FALSE)
      }
      results_list[[length(results_list) + 1]] <- cell_results
    }, error = function(e) {
      # `<<-` because this assignment happens inside the error handler's own
      # scope but must update the loop-level `skipped` list.
      skipped[[length(skipped) + 1]] <<- tibble(
        study = cell$study, contrast = cell$contrast, direction = cell$direction,
        n_query = length(cell$query), reason = conditionMessage(e))
    })
  }

  combined <- if (length(results_list)) bind_rows(results_list) else tibble()
  skipped_df <- if (length(skipped)) bind_rows(skipped) else
    tibble(study = character(), contrast = character(), direction = character(),
           n_query = integer(), reason = character())

  if (isTRUE(write)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(combined, file.path(out_dir, paste0(outcome_group, "_", arm_set, "_all_cells.csv")),
              row.names = FALSE)
    write.csv(skipped_df, file.path(out_dir, "skipped_cells.csv"), row.names = FALSE)
  }
  list(results = combined, skipped = skipped_df, out_dir = out_dir)
}
