# =============================================================================
# fig6C-proteomics.R
#
# Builds Fig 6C, the GO Biological-Process over-representation of the milk proteome
# as a volcano in the style of Panels A and B: x = signed fold enrichment,
# y = -log10(P), colour and shape = study when P < 0.05 (grey otherwise), P < 0.05
# and Q < 0.05 threshold lines, and the top 3 FDR-significant terms per study x
# direction labelled. The GO run counts each UniProt protein once; under it only
# down-regulated terms reach FDR significance (MISAME-III telomere/DNA, Mumta-LW
# nucleotide/NAD/energy terms).
#
# Inputs:  results/proteomics_go_uniprot.csv (src/2 analysis/55-proteomics-go-uniprot.R;
#            also the source of Table S7)
# Outputs: figures/figure6_panelC_proteomics_go.png
# =============================================================================
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})
root <- here::here()
source(file.path(root, "figure-scripts/manuscript_figures/study_colors.R"))  # shared study colours and shapes
source(file.path(root, "figure-scripts/0_figure-functions.R"))               # theme_imic(), imic_logp_title
source(file.path(root, "figure-scripts/manuscript_figures/fig6_layout.R"))    # FIG6_* geometry + print-scale text sizes

d <- fread(file.path(root, "results/proteomics_go_uniprot.csv"))
d[, `:=`(fe = signed_fold_enrichment, y = -log10(pvalue), sig = p.adjust < 0.05)]
# Only MISAME-III and Mumta-LW have proteomics; index the shared study map by name so
# each keeps its colour (orange / purple), never by order.
studycols <- imic_study_cols[c("MISAME-III", "Mumta-LW")]
# Same tiers as Panel A and Figs 3B/5B: every P < 0.05 term filled in its study colour
# and symbol, the rest grey. The green line marks the FDR frontier.
d[, study_col := fifelse(pvalue < 0.05, study, "Not Significant")]
setorder(d[, draw := fifelse(study_col == "Not Significant", 0L, 1L)], draw)  # coloured on top
q_line <- suppressWarnings(min(d[sig == TRUE, y]))     # -log10(P) at the BH-FDR 0.05 frontier
# GO-term abbreviations, as in the original R Markdown proteomics analysis.
abbr_proteomics <- function(x) {
  reps <- c(
    "nicotinamide nucleotide metabolism" = "Nicotinamide nt metab",
    "pyridine-containing compound metabolism" = "Pyridine compound metab",
    "purine ribonucleoside diphosphate catabolic process" = "Purine ribonucleoside DP catab",
    "purine nucleoside diphosphate catabolic process" = "Purine nucleoside DP catab",
    "pyridine nucleotide catabolic process" = "Pyridine nt catab",
    "nucleotide catabolic process" = "Nucleotide catab",
    "ADP catabolic process" = "ADP catab", "glycolytic process" = "Glycolysis",
    "glucose metabolism" = "Glucose metab", "nucleotide metabolism" = "Nucleotide metab",
    "positive regulation of protein modification process" = "Pos reg protein modification",
    "negative regulation" = "Neg reg", "positive regulation" = "Pos reg",
    "ribonucleoprotein complex" = "RNP complex", "protein localization" = "Protein loc",
    "intracellular" = "Intracell", "biosynthetic process" = "Biosynthesis",
    "catabolic process" = "Catabolism", "metabolic process" = "Metab",
    "metabolism" = "Metab", "regulation" = "Reg", "organization" = "Org")
  for (i in seq_along(reps)) x <- gsub(names(reps)[i], reps[[i]], x, ignore.case = TRUE)
  x <- tools::toTitleCase(trimws(gsub("\\s+", " ", x)))
  x <- gsub("\\bAdp\\b", "ADP", x); x <- gsub("\\bDp\\b", "DP", x); gsub("\\bRnp\\b", "RNP", x)
}
d[, lab := abbr_proteomics(description_sentence)]
# Label FDR-significant terms: top 3 per study x direction by |fold enrichment| (a
# half-column panel has room for few labels).
labs <- d[sig == TRUE][order(-abs(fe))][, .SD[!duplicated(lab)], by = .(study, direction)][
  , head(.SD, 3), by = .(study, direction)]
p_cap_x <- fig6_caption_x(d$fe, d$y, -log10(0.05))
q_cap_x <- fig6_caption_x(d$fe, d$y, q_line)
cap_obst <- rbind(fig6_caption_obstacles(p_cap_x, -log10(0.05), d$fe),
                  fig6_caption_obstacles(q_cap_x, q_line, d$fe))   # empty labels: repel only
lab_rep <- rbind(labs, data.table(cap_obst, lab = "", study = "Not Significant"), fill = TRUE)

# FIG6_Y_MAX / FIG6_PANEL_*: fig6_layout.R. No legend of its own: one shared Fig 6 A-C
# key (Panel A's) sits under this panel, and A/B/C share one panel size.
pal <- c(studycols, "Not Significant" = "grey75")
tier_shapes <- imic_shapes_for(names(pal))
p <- ggplot(d, aes(fe, y)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
  # line colours/weights and captions as in Figs 3B and 5B
  geom_hline(yintercept = -log10(0.05), colour = "#BAB0AC", linewidth = 0.4) +
  geom_hline(yintercept = q_line, colour = "#59A14F", linewidth = 0.4) +
  geom_point(aes(colour = study_col, shape = study_col), size = 2, alpha = 0.85) +  # study symbol = colourblind cue
  # captions on the emptiest stretch of each line (fig6_caption_x(), fig6_layout.R), on a
  # borderless white box above the points: if no stretch of a line is empty, the box
  # hides the few points under the caption
  annotate("label", x = p_cap_x, y = -log10(0.05), label = 'italic(P)*"-value" < 0.05', parse = TRUE,
           hjust = 0, vjust = -0.3, size = 2.4, colour = "#8A8580", fill = "white",
           border.colour = NA, label.padding = unit(0.4, "mm")) +
  annotate("label", x = q_cap_x, y = q_line, label = 'italic(Q)*"-value" < 0.05', parse = TRUE,
           hjust = 0, vjust = -0.3, size = 2.4, colour = "#3C8C3C", fill = "white",
           border.colour = NA, label.padding = unit(0.4, "mm")) +
  # size 2.2 / label.padding 0.1 (and base_size 8 below), as in fig6B-mummichog.R
  geom_label_repel(data = lab_rep, aes(label = lab, colour = study), size = 2.2,
                  max.overlaps = Inf, min.segment.length = 0, box.padding = 0.3,
                  label.padding = 0.1, label.size = 0.15, fill = "white",  # white box (matches Fig 6A)
                  seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = pal,
                      breaks = c("MISAME-III", "Mumta-LW", "Not Significant"),
                      name = "Study") +
  scale_shape_manual(values = tier_shapes,
                     breaks = c("MISAME-III", "Mumta-LW", "Not Significant"),
                     name = "Study") +
  scale_y_continuous(limits = c(0, FIG6_Y_MAX), breaks = seq(0, FIG6_Y_MAX, 1)) +  # shared Fig 6 A/B/C height
  labs(x = "Fold Enrichment (signed by direction)",
       y = imic_logp_title) +
  theme_imic(base_size = 8) +   # shared theme (Helvetica, font-size floors); matches Fig 3A/5A
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "none", legend.key.size = unit(0.35, "cm"),
        # print at Figs 3/5 sizes; .x/.y because theme_imic() sizes both children explicitly
        axis.text.x = element_text(size = fig6_pt(FIG6_TICK_PT, FIG6_ABC_SCALE)),
        axis.text.y = element_text(size = fig6_pt(FIG6_TICK_PT, FIG6_ABC_SCALE)),
        axis.title = element_text(size = fig6_pt(FIG6_TITLE_PT, FIG6_ABC_SCALE)),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        plot.margin = margin(4, 6, 4, 4))

# Left-column half-page sub-panel of the full-page (7.25 in) Fig 6 (A/B/C left, D right).
ggsave(file.path(root, "figures/figure6_panelC_proteomics_go.png"),
       p, width = FIG6_PANEL_W_MM, height = FIG6_PANEL_H_MM, units = "mm", dpi = 300, bg = "white",  # 105 mm wide, shared A/B/C height
       device = ragg::agg_png)
cat("wrote figures/figure6_panelC_proteomics_go.png |", nrow(d), "terms,",
    sum(d$sig), "FDR-sig,", nrow(labs), "labelled | up FDR-sig:",
    d[sig == TRUE & direction == "Upregulated", .N], "\n")
