# =============================================================================
# fig3B-pathway.R
#
# Builds Fig 3B: KEGG pathway impact (x) against enrichment -log10(P) (y) for the
# primary B-vitamin outcomes, one point per pathway x study x visit, coloured and
# shaped by study (grey when P >= 0.05). Pathways nominally significant in at least
# two study x visit cells are grouped by an ellipse and labelled. The pathway run
# uses the KEGG library of the MetaboAnalyst web tool (metpa, pathlib = "kegg"),
# which reproduces the submitted impact values; raw P-values can differ slightly
# from the web exports because of compound-name database drift. The PNG is embedded
# in Fig 3 by fig3-primary-volcano-composite.R.
#
# Inputs:  results/metaboanalyst/primary_pathway_local/primary_pathway_all_cells.csv
#            (src/metaboanalyst/run-primary-pathway-local.R; also Table S2)
# Outputs: figures/figure3_panelB_msea.png
# Run from the repo root (paths are relative).
# =============================================================================
suppressMessages({ library(dplyr); library(ggplot2); library(ggrepel); library(stringr) })
source("figure-scripts/manuscript_figures/study_colors.R")   # shared study colours and shapes
source("figure-scripts/0_figure-functions.R")                # imic_logp_title

LOCAL_CSV <- "results/metaboanalyst/primary_pathway_local/primary_pathway_all_cells.csv"
OUT_EMBED <- "figures/figure3_panelB_msea.png"

# Fig 3 embeds this PNG (saved 105 x 99 mm) height-limited in a 0.32 x 10.63 in cell,
# i.e. at ~87% of its saved size. Theme text sizes are divided by this so Panel B's
# axis text and titles print at Panel A's sizes (8 pt) in the composite.
EMBED_SCALE <- (0.32 * 10.63) / (297 / 3 / 25.4)

# Panel style: solid P < 0.05 (grey) and Q < 0.05 (green) reference lines with
# left-anchored captions, white boxed labels and ellipses, integer y ticks, x breaks
# at 0.25, a bottom-left legend, and a redundant study shape (circle / triangle /
# square) alongside the study colour.
.line_caption_obstacles <- function(y_fdr_line, caption_x, caption_hjust) {
  ys <- c(-log10(0.05), if (is.finite(y_fdr_line)) y_fdr_line) + 0.18
  x0 <- caption_x - caption_hjust * 0.2   # captions are ~0.2 impact units wide
  expand.grid(impact = seq(x0, x0 + 0.2, by = 0.025), logP = ys) %>%
    mutate(lab = "", point_color = "Not Significant")
}

# caption_x / caption_hjust place the "P-value < 0.05" / "Q-value < 0.05" line captions
# in an empty stretch of the panel.
.plot_3b <- function(tab, ellipse_df, label_df, cols, fdr_p_thr, y_fdr_line,
                     caption_x = 0.045, caption_hjust = 0) {
  ggplot(tab, aes(impact, logP)) +
    geom_hline(yintercept = -log10(0.05), color = "#BAB0AC", linewidth = 0.4) +   # solid P<0.05
    { if (is.finite(y_fdr_line)) geom_hline(yintercept = y_fdr_line, color = "#59A14F", linewidth = 0.4) } +  # solid Q<0.05
    annotate("text", x = caption_x, y = -log10(0.05) + 0.18, hjust = caption_hjust, size = 2.3, color = "#8A8580",
             parse = TRUE, label = "italic(P)*\"-value\" < 0.05") +
    # green line = FDR threshold, labelled by what it means (Q < 0.05) rather than the
    # raw P it falls at
    { if (is.finite(y_fdr_line)) annotate("text", x = caption_x, y = y_fdr_line + 0.18, hjust = caption_hjust, size = 2.3,
             color = "#3C8C3C", parse = TRUE,
             label = "italic(Q)*\"-value\" < 0.05") } +
    { if (nrow(ellipse_df) >= 2) ggforce::geom_mark_ellipse(data = ellipse_df, aes(group = pathway),
             expand = unit(2, "mm"), colour = "black", linewidth = 0.3) } +
    geom_point(aes(color = point_color, shape = point_color), size = 1.6, alpha = 0.9) +
    # unlabelled points along the two reference-line captions: ggrepel keeps the
    # pathway labels clear of every point in its data, including ones labelled ""
    geom_label_repel(data = bind_rows(label_df, .line_caption_obstacles(y_fdr_line, caption_x, caption_hjust)),
             aes(label = lab, color = point_color),
             size = 2.2, fill = "white", label.size = 0.2, box.padding = 0.2,
             point.padding = 0.3, segment.size = 0.35, max.overlaps = Inf,
             show.legend = FALSE, min.segment.length = 0, seed = 123,
             # the default 0.5 s time limit left labels overlapping; max.time = Inf stops
             # on the iteration count alone, so placement cannot depend on machine speed
             max.time = Inf, max.iter = 1e5, force = 3) +
    scale_color_manual(values = cols, name = "Study") +
    scale_shape_manual(values = imic_shapes_for(names(cols)), name = "Study") +   # redundant CVD cue
    scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.05, 1.05)) +
    scale_y_continuous(breaks = seq(0, ceiling(max(tab$logP, na.rm = TRUE)), 1)) +
    labs(x = "Pathway Impact", y = imic_logp_title) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    # Panel A's look (theme_imic): no y ticks and no horizontal gridlines, vertical
    # gridlines only. Explicit sizes, divided by EMBED_SCALE so they print at Panel
    # A's 8 pt axis text/titles in the composite.
    theme(panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks.y       = element_blank(),
          # thin black border, as on Panel A and the other panels
          panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.3 / EMBED_SCALE),
          axis.text        = element_text(size = 8 / EMBED_SCALE),
          axis.title       = element_text(size = 8 / EMBED_SCALE),
          legend.text      = element_text(size = 7 / EMBED_SCALE),
          legend.title     = element_text(size = 8 / EMBED_SCALE),
          legend.position = "bottom", legend.justification = "left", legend.box = "horizontal",
          # anchor the 7 pt legend to the plot's left edge (not the panel's) and
          # tighten the key spacing so "Not Significant" fits
          legend.location = "plot", legend.key.spacing.x = unit(1.5, "mm"),
          legend.margin = margin(0, 0, 0, 0)) +
    guides(color = guide_legend(override.aes = list(size = 2.5)))   # shape legend merges into this one
}

build <- function() {
  # local KEGG pathway run: already tidy (study, tp, pathway, impact, raw_p, fdr).
  tab <- read.csv(LOCAL_CSV, check.names = FALSE, stringsAsFactors = FALSE) %>%
    transmute(study, tp, pathway, impact = as.numeric(impact),
              raw_p = as.numeric(raw_p), fdr = as.numeric(fdr))

  tab <- tab %>% mutate(
    logP = -log10(raw_p),
    is_sig = raw_p < 0.05,
    point_color = ifelse(is_sig, study, "Not Significant"),
    short = pathway %>%
      str_replace(regex("Nicotinate and nicotinamide metabolism", ignore_case = TRUE), "NAD/NAM Met") %>%
      str_replace(regex("Riboflavin metabolism", ignore_case = TRUE), "Riboflavin Met") %>%
      str_replace(regex("Thiamine metabolism", ignore_case = TRUE), "Thiamine Met") %>%
      str_replace(regex("Vitamin B6 metabolism", ignore_case = TRUE), "Vit B6 Met") %>%
      str_replace(regex("Pantothenate and CoA biosynthesis", ignore_case = TRUE), "Pantothenate Met") %>%
      str_replace("metabolism", "Met"),
    lab = ifelse(is_sig, paste0(short, " (", tp, ")"), NA))

  cols <- c(imic_study_cols, "Not Significant"="grey75")
  rep_df <- tab %>% filter(is_sig) %>% add_count(pathway, name = "pc") %>% filter(pc >= 2)

  # Green FDR reference line = the highest raw p still FDR-significant (pooled across
  # cells), matching the Q line on the MSEA panels (render_msea_panelB.R). FDR is the
  # per-cell Benjamini-Hochberg value from the local pathway run.
  fdr_p_thr  <- suppressWarnings(max(tab$raw_p[tab$fdr < 0.05], na.rm = TRUE))
  y_fdr_line <- if (is.finite(fdr_p_thr)) -log10(fdr_p_thr) else NA_real_

  # combined-arm panel; labels the pathways significant in >= 2 study x time cells.
  p <- .plot_3b(tab, ellipse_df = rep_df, label_df = rep_df, cols = cols,
                fdr_p_thr = fdr_p_thr, y_fdr_line = y_fdr_line)

  # Panel size 105 x 99 mm (210/2 x 297/3, near-square), as submitted; Fig 5C uses
  # the same size.
  ggsave(OUT_EMBED, p, width = 210/2, height = 297/3, units = "mm", dpi = 300, device = ragg::agg_png)
  cat("wrote 3B from LOCAL KEGG pathway run | cells:", n_distinct(paste(tab$study, tab$tp)),
      "| impact range:", round(range(tab$impact, na.rm = TRUE), 3), "| repeated-pathway groups:",
      n_distinct(rep_df$pathway), "\n")
  invisible(p)
}

if (sys.nframe() == 0) invisible(build())
