# =============================================================================
# figS10-compartment-tracking.R
#
# Supplementary Figure S10: cross-compartment tracking of BEP-associated
# UP-REGULATED metabolite features and their BEP-supplement origin (MISAME-III).
#
# WHY THIS SCRIPT EXISTS
#   The shipped Fig. S10 (figures/cross_compartment/trenton_xcompartment_tracking.png,
#   2026-08-12) is an orphan binary: it came from Trenton's analysis and no script in
#   this repo produced it, so it could not be rebuilt when its inputs changed.
#
# METHOD FIDELITY -- READ THIS BEFORE CHANGING THE GROUPING
#   This script reproduces the grouping rule in Trenton's
#   `src/trenton-manuscript-analyses/Compartment Tracking.Rmd`, VERBATIM:
#
#     lines 1743-1749 : tracking_name = coalesce(putative_name, paste0("m/z ", round(mz, 4)))
#     lines 2556-2562 : group_by(direction, tracking_name) %>%
#                         filter(n_distinct(compartment_timepoint) >= 2)
#
#   Two consequences, both deliberate and both easy to get wrong:
#
#   (1) The unit of analysis is the FEATURE'S OWN putative name, not a connected
#       component of matched features. Two features linked only by 25-ppm mass are
#       plotted as SEPARATE rows if their putative annotations differ. This is why
#       D-lactate (m/z 135.03) and N1-methyl-2-pyridone-5-carboxamide (m/z 153.066)
#       appear under their own names rather than being absorbed by their stronger
#       isobars (threonic acid, N(pi)-methyl-L-histidine).
#
#       `src/metaboanalyst/run-compartment-pathway.R` builds a union-find
#       `tracking_group` instead and names each component after its lowest-P member.
#       That is a GENERALISATION of his method, not his method: it merges the isobars
#       above and drops those two compounds from the figure and from the Table S11
#       compound list. `results/compartment_tracking/linked_upregulated_plot_data.csv`
#       is that union-find output and is therefore NOT used here.
#
#   (2) The inclusion threshold is >= 2 distinct COMPARTMENT x TIMEPOINT CELLS, which
#       a feature can satisfy inside ONE compartment. That is why the shipped figure's
#       top rows (m/z 286.2021, m/z 286.1984, m/z 218.103) each show a single
#       compartment. The printed caption's "linked across at least two of maternal
#       blood, milk, and infant blood" describes a SUBSET of the rows, not the
#       inclusion rule. The figure keeps his rule; the compartment-spanning counts the
#       caption needs are written to results/tables/figS10_tracking_counts.csv.
#
#   Set GROUPING = "union-find" to render the alternative for comparison.
#
# Inputs : results/supplement_status_fdr_features.csv  (script 54)
#          results/compartment_tracking/metabolite_pathways_local.csv  (pathway labels)
#          src/trenton-ports/compartment-tracking.R    (his pairwise matcher, faithful)
# Outputs: figures/cross_compartment/figureS10_compartment_tracking.png
#          results/tables/figS10_tracking_counts.csv
#          results/compartment_tracking/trenton_linked_upregulated.csv  (row-level data)
#          results/compartment_tracking/trenton_pathway_compound_list.csv (Table S11 query)
#
# Run from repo root:
#   Rscript figure-scripts/manuscript_figures/figS10-compartment-tracking.R
# =============================================================================
suppressMessages({ library(dplyr); library(tidyr); library(stringr); library(ggplot2)
                   library(forcats); library(readr) })
source("figure-scripts/0_figure-functions.R")
source("src/trenton-ports/compartment-tracking.R")   # ct_compartment_pairs(), ct_normalized_name()

GROUPING <- "trenton"          # "trenton" (his name-based rule) | "union-find"
SIGFEAT  <- "results/supplement_status_fdr_features.csv"
# pathway labels: the "trenton" grouping's own KEGG run (written by
# run-compartment-pathway-trenton.R, run AFTER this script); the union-find
# alt keeps reading the union-find run's file (run-compartment-pathway.R) so
# the labels always match whichever grouping produced the row's tracking_name.
PATHS_TRENTON    <- "results/compartment_tracking/metabolite_pathways_trenton.csv"
PATHS_UNIONFIND  <- "results/compartment_tracking/metabolite_pathways_local.csv"
OUT      <- "figures/cross_compartment/figureS10_compartment_tracking.png"
COUNTS   <- "results/tables/figS10_tracking_counts.csv"
ROWDATA  <- "results/compartment_tracking/trenton_linked_upregulated.csv"
CMPDLIST <- "results/compartment_tracking/trenton_pathway_compound_list.csv"

COMP_LEVELS <- c("Maternal plasma", "Maternal VAMS", "Milk", "Infant VAMS")
COMP_COLS   <- c("Maternal plasma" = "#4E79A7", "Maternal VAMS" = "#59A14F",
                 "Milk" = "#F28E2B", "Infant VAMS" = "#E15759")
# his compartment x timepoint ordering (Compartment Tracking.Rmd lines 1705-1740)
CT_LEVELS <- c("Maternal plasma · incl","Maternal plasma · tri3","Maternal plasma · pn12",
               "Maternal VAMS · tri3","Maternal VAMS · pn56","Milk · 1421d","Milk · pn12",
               "Milk · pn34","Infant VAMS · acco","Infant VAMS · pn12","Infant VAMS · pn34",
               "Infant VAMS · pn56")
# maternal plasma + maternal VAMS are two assays of ONE compartment; the printed
# caption's three compartments are maternal blood / milk / infant blood.
POOL3 <- c("Maternal plasma" = "Maternal blood", "Maternal VAMS" = "Maternal blood",
           "Milk" = "Milk", "Infant VAMS" = "Infant blood")
SUPP_LABS   <- c(`TRUE` = "Detected in supplement", `FALSE` = "Not reliably detected in supplement")
SUPP_SHAPES <- c("Detected in supplement" = 16, "Not reliably detected in supplement" = 17)
# junk annotation carried by the upstream name map; excluded from the pathway query
JUNK_NAMES <- c("Acid")

is_true <- function(x) { if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes") }

build_figS10 <- function(grouping = GROUPING, out_png = OUT) {
  sig <- read_csv(SIGFEAT, show_col_types = FALSE) %>% filter(!is.na(mz), !is.na(effect_size))

  # --- his pairwise linkage (ported verbatim in src/trenton-ports) -------------
  pairs   <- ct_compartment_pairs(sig)
  tracked <- pairs %>% filter(same_direction) %>%
    select(tracking_id_1, tracking_id_2) %>%
    pivot_longer(everything(), values_to = "tracking_id") %>% distinct(tracking_id)

  # shared grouping rule (src/trenton-ports/compartment-tracking.R), computed
  # once over the full significant set -- "trenton" = his name-based
  # tracking_name; "union-find" = the connected-component alt, kept for
  # comparison only.
  sig_tid <- sig %>% mutate(tracking_id = paste(compartment, timepoint, feature, sep = "__"))
  group_lookup <- if (grouping == "trenton") {
    ct_trenton_track(sig_tid) %>% select(tracking_id, tracking_name)
  } else {
    ct_tracking_components(sig) %>% select(tracking_id, tracking_name = tracking_group)
  }

  d <- sig_tid %>%
    semi_join(tracked, by = "tracking_id") %>%
    left_join(group_lookup, by = "tracking_id") %>%
    mutate(compartment_timepoint = factor(paste(compartment, timepoint, sep = " · "),
                                          levels = CT_LEVELS),
           comp3     = unname(POOL3[compartment]),
           reliably  = is_true(supplement_reliably_detected))

  # --- his ">= 2 positions" rule, within direction -----------------------------
  up <- d %>% filter(direction == "Upregulated", !is.na(tracking_name), tracking_name != "") %>%
    group_by(tracking_name) %>%
    filter(n_distinct(compartment_timepoint) >= 2) %>%
    arrange(compartment_timepoint, raw_p, .by_group = TRUE) %>%
    distinct(tracking_name, compartment_timepoint, .keep_all = TRUE) %>%
    ungroup()

  per <- up %>% group_by(tracking_name) %>%
    summarise(n_cells = n_distinct(compartment_timepoint),
              n_comp3 = n_distinct(comp3), reliably = any(reliably), .groups = "drop")
  span2 <- per %>% filter(n_comp3 >= 2)

  dir.create(dirname(COUNTS), recursive = TRUE, showWarnings = FALSE)
  write_csv(tibble(
    metric = c("grouping_rule", "n_rows_in_figure", "n_features_ge2_compartments",
               "n_reliably_supplement_matched", "pct_reliably_supplement_matched",
               "n_features_all3_compartments", "features_all3_compartments"),
    value  = c(grouping, as.character(nrow(per)), as.character(nrow(span2)),
               as.character(sum(span2$reliably)),
               sprintf("%.1f", 100 * mean(span2$reliably)),
               as.character(sum(per$n_comp3 >= 3)),
               paste(sort(per$tracking_name[per$n_comp3 >= 3]), collapse = "; "))), COUNTS)
  write_csv(up, ROWDATA)
  write_csv(tibble(compound = sort(per$tracking_name[!grepl("^m/z ", per$tracking_name) &
                                                     !per$tracking_name %in% JUNK_NAMES])),
            CMPDLIST)

  # --- pathway-annotated y-axis labels ----------------------------------------
  # Two-pass build: this script must run FIRST to emit CMPDLIST (the pathway
  # query list), then src/metaboanalyst/run-compartment-pathway-trenton.R runs
  # the KEGG ORA and writes PATHS_TRENTON, then this script runs a SECOND time
  # to pick up the labels for the final image. On the first pass PATHS_TRENTON
  # does not exist yet -- the figure still renders, just without pathway text.
  paths_file <- if (grouping == "trenton") PATHS_TRENTON else PATHS_UNIONFIND
  MAX_PATHS <- 2L
  paths <- if (file.exists(paths_file)) read_csv(paths_file, show_col_types = FALSE) else
    tibble(tracking_name = character(), pathway = character(), raw_p = numeric())
  paths <- paths %>% arrange(raw_p) %>%
    mutate(line = paste0(pathway, " (P = ",
                         ifelse(raw_p < 0.001, "<0.001", sprintf("%.3f", raw_p)), ")")) %>%
    distinct(tracking_name, line, .keep_all = TRUE) %>%
    group_by(tracking_name) %>% slice_head(n = MAX_PATHS) %>% ungroup()
  axis_labels <- paths %>% group_by(tracking_name) %>%
    summarise(pathway_text = paste(line, collapse = "\n"), .groups = "drop")

  p_d <- up %>% left_join(axis_labels, by = "tracking_name") %>%
    mutate(row_label = ifelse(is.na(pathway_text), as.character(tracking_name),
                              paste0(tracking_name, "\n", pathway_text)),
           row_label = fct_reorder(row_label, effect_size, .fun = max),
           compartment = factor(compartment, levels = COMP_LEVELS),
           supp = factor(unname(SUPP_LABS[as.character(reliably)]), levels = names(SUPP_SHAPES)))
  seg <- p_d %>% group_by(row_label) %>%
    summarise(xmin = min(effect_size), xmax = max(effect_size), .groups = "drop")

  p <- ggplot(p_d, aes(x = effect_size, y = row_label)) +
    geom_segment(data = seg, aes(x = xmin, xend = xmax, y = row_label, yend = row_label),
                 inherit.aes = FALSE, colour = "grey65", linewidth = 0.7) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_point(aes(colour = compartment, shape = supp), size = 2.6, alpha = 0.95) +
    scale_colour_manual(values = COMP_COLS, drop = FALSE, name = "Compartment") +
    scale_shape_manual(values = SUPP_SHAPES, drop = FALSE, name = "BEP supplement") +
    labs(x = "Average Treatment Effect", y = NULL) +
    guides(shape = guide_legend(order = 1, nrow = 2, override.aes = list(size = 2.4)),
           colour = guide_legend(order = 2, nrow = 2, override.aes = list(size = 2.4))) +
    theme_imic(base_size = 9) +
    theme(axis.text.y = element_text(size = 6.5, lineheight = 0.9, hjust = 1),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          axis.line.x = element_line(colour = "black", linewidth = 0.3),
          legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
          legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
          legend.key.size = unit(3.5, "mm"), legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "mm"), plot.margin = margin(6, 10, 6, 6))

  n_lines <- nlevels(p_d$row_label) + nrow(paths %>% filter(tracking_name %in% up$tracking_name))
  h <- min(9.5, 1.7 + 0.165 * n_lines)
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_png, p, width = 5.0, height = h, units = "in", dpi = 300,
         bg = "white", device = ragg::agg_png)

  cat("wrote", out_png, "  [grouping:", grouping, "]\n")
  cat("  rows in figure (>=2 compartment x timepoint cells):", nrow(per), "\n")
  cat("  spanning >=2 of {maternal blood, milk, infant blood}:", nrow(span2), "\n")
  cat("  ... of which reliably supplement-matched:", sum(span2$reliably),
      sprintf("(%.1f%%)\n", 100 * mean(span2$reliably)))
  cat("  spanning all three:", sum(per$n_comp3 >= 3), "->",
      paste(sort(per$tracking_name[per$n_comp3 >= 3]), collapse = ", "), "\n")
  cat("  pathway query compounds ->", CMPDLIST, "\n")
  invisible(p)
}

if (sys.nframe() == 0) invisible(build_figS10())
