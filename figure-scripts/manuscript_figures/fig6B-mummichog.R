# =============================================================================
# fig6B-mummichog.R
#
# Builds Fig 6B, the mummichog pathway plot for the untargeted milk metabolome:
# signed pathway coverage (x) against -log10(P) (y), one point per pathway x study x
# visit x ionization mode, coloured and shaped by study when P < 0.05 (grey
# otherwise), with P < 0.05 (grey) and Q < 0.05 (green) lines and a vertical line
# at 0. The pathways discussed in the Results are labelled with a "+"/"-"
# ionization suffix and the visit.
#
# Inputs:  results/metaboanalyst/mummichog/milk_mummichog_pathways.csv
#            (src/metaboanalyst/run-milk-mummichog.R; Table S6)
# Outputs: figures/figure6_panelB_mummichog.png
# =============================================================================
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})
root <- here::here()
source(file.path(root, "figure-scripts/manuscript_figures/study_colors.R"))  # shared study colours and shapes
source(file.path(root, "figure-scripts/0_figure-functions.R"))               # theme_imic(), imic_logp_title
source(file.path(root, "figure-scripts/manuscript_figures/fig6_layout.R"))    # FIG6_* geometry + print-scale text sizes

d <- fread(file.path(root, "results/metaboanalyst/mummichog/milk_mummichog_pathways.csv"))
setnames(d, c("P-value", "Time Point", "Ionization Mode"), c("p_value", "timepoint", "ionization"))
d <- d[!is.na(`Enrichment Ratio`) & !is.na(p_value)]
d[, `:=`(
  fe       = `Enrichment Ratio`,                       # already signed by direction (neg = down)
  y        = -log10(p_value),
  study    = fifelse(Study == "MISAME", "MISAME-III", Study),
  ion_sign = fifelse(ionization == "Positive", "+", "\u2013"),   # en dash: easier to see than "-"
  sig      = FDR < 0.05,
  # nominal (before-FDR) significance: raw P < 0.05 but not FDR-significant.
  nom_sig  = (p_value < 0.05) & !(FDR < 0.05)
)]
# Shared study colours (ELICIT blue / MISAME orange / Mumta purple); index by name.
studycols <- c(imic_study_cols, "Not Significant" = "grey75")
# Same tiers as Panel A and Figs 3B/5B: every P < 0.05 pathway filled in its study colour
# and symbol, the rest grey. The green line marks the FDR frontier.
d[, study_col := fifelse(p_value < 0.05, study, "Not Significant")]
setorder(d[, draw := fifelse(study_col == "Not Significant", 0L, 1L)], draw)  # coloured on top
tier_shapes <- imic_shapes_for(names(studycols))
d[, lab := fifelse(sig, paste0(Pathway, ion_sign, " (", timepoint, ")"), NA_character_)]
# Label the pathways discussed in the Results (each is FDR-significant here), one label
# per distinct pathway (the most extreme |fe| instance). "Prostaglandin formation from
# arachidonate" is FDR-significant only in the up direction (fe > 0), outside the
# coordinated reduction the text describes, so it is not labelled.
discussed <- c("Arachidonic acid metabolism", "Linoleate metabolism",
               "De novo fatty acid biosynthesis", "Fatty acid activation",
               "Fatty Acid Metabolism", "Leukotriene metabolism")
lab_d <- d[sig == TRUE & Pathway %in% discussed][order(-abs(fe))][!duplicated(Pathway)]
q_line <- suppressWarnings(min(d[sig == TRUE, y]))
p_cap_x <- fig6_caption_x(d$fe, d$y, -log10(0.05))
q_cap_x <- fig6_caption_x(d$fe, d$y, q_line)
cap_obst <- rbind(fig6_caption_obstacles(p_cap_x, -log10(0.05), d$fe),
                  fig6_caption_obstacles(q_cap_x, q_line, d$fe))   # empty labels: repel only
lab_rep <- rbind(lab_d, data.table(cap_obst, lab = "", study_col = "Not Significant"), fill = TRUE)
# FIG6_Y_MAX / FIG6_PANEL_*: fig6_layout.R. No legend of its own: one shared Fig 6 A-C
# key (Panel A's) sits under Panel C, and A/B/C share one panel size.

p <- ggplot(d, aes(fe, y, colour = study_col)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
  # line colours/weights and captions as in Figs 3B and 5B
  geom_hline(yintercept = -log10(0.05), colour = "#BAB0AC", linewidth = 0.4) +
  geom_hline(yintercept = q_line, colour = "#59A14F", linewidth = 0.4) +
  geom_point(aes(shape = study_col), size = 2, alpha = 0.85) +   # study symbol = colourblind cue
  # captions on the emptiest stretch of each line (fig6_caption_x(), fig6_layout.R), on a
  # borderless white box above the points: in Panel B points sit just above the P line
  # at every x, so no stretch is empty; the box hides the few points under the caption
  annotate("label", x = p_cap_x, y = -log10(0.05), label = 'italic(P)*"-value" < 0.05', parse = TRUE,
           hjust = 0, vjust = -0.3, size = 2.4, colour = "#8A8580", fill = "white",
           border.colour = NA, label.padding = unit(0.4, "mm")) +
  annotate("label", x = q_cap_x, y = q_line, label = 'italic(Q)*"-value" < 0.05', parse = TRUE,
           hjust = 0, vjust = -0.3, size = 2.4, colour = "#3C8C3C", fill = "white",
           border.colour = NA, label.padding = unit(0.4, "mm")) +
  # size 2.2 / label.padding 0.1 (and base_size 8 below) to match the smaller, boxless
  # repel-label text of Figs 3A/5A
  geom_label_repel(data = lab_rep, aes(label = lab), size = 2.2, na.rm = TRUE,
                  max.overlaps = Inf, min.segment.length = 0, box.padding = 0.5,
                  label.padding = 0.1, label.size = 0.15, fill = "white",  # white box (matches Fig 6A)
                  force = 6, nudge_y = 0.6, segment.size = 0.25,
                  seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = studycols,
                      breaks = c("ELICIT", "MISAME-III", "Mumta-LW", "Not Significant"),
                      name = "Study") +
  scale_shape_manual(values = tier_shapes,
                     breaks = c("ELICIT", "MISAME-III", "Mumta-LW", "Not Significant"),
                     name = "Study") +
  scale_y_continuous(limits = c(0, FIG6_Y_MAX), breaks = seq(0, FIG6_Y_MAX, 1)) +
  # x = Table S6's "Enrichment Ratio", which is overlap size / pathway size (the share of
  # the pathway's compounds among the significant features; run-milk-mummichog.R),
  # signed by direction. It is not a fold enrichment (hits / expected) like Panels A and C;
  # "pathway coverage" is the Methods' name for it.
  labs(x = "Pathway Coverage (signed by direction)", y = imic_logp_title) +
  theme_imic(base_size = 8) +   # shared theme (Helvetica, font-size floors); matches Fig 3A/5A
  guides(colour = guide_legend(nrow = 2, byrow = TRUE),
         shape  = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "none", legend.key.size = unit(0.35, "cm"),
        # print at Figs 3/5 sizes; .x/.y because theme_imic() sizes both children explicitly
        axis.text.x = element_text(size = fig6_pt(FIG6_TICK_PT, FIG6_ABC_SCALE)),
        axis.text.y = element_text(size = fig6_pt(FIG6_TICK_PT, FIG6_ABC_SCALE)),
        axis.title = element_text(size = fig6_pt(FIG6_TITLE_PT, FIG6_ABC_SCALE)),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        plot.margin = margin(3, 6, 3, 3))

# Left-column half-page sub-panel of the full-page (7.25 in) Fig 6 (A/B/C left, D right).
ggsave(file.path(root, "figures/figure6_panelB_mummichog.png"),
       p, width = FIG6_PANEL_W_MM, height = FIG6_PANEL_H_MM, units = "mm", dpi = 300, bg = "white")  # 105 mm wide, shared A/B/C height
cat("wrote figures/figure6_panelB_mummichog.png |", nrow(d), "pathways,",
    sum(d$sig), "FDR-sig,", nrow(lab_d), "labelled\n")
