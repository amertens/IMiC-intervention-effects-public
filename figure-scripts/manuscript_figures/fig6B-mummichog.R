# fig6B-mummichog.R
# Fig 6 Panel B: Mummichog pathway plot for the untargeted milk metabolome,
# RESTYLED to match the SUBMITTED figure: fold-enrichment CIRCLES coloured by
# study (canonical colours: ELICIT blue, MISAME-III orange, Mumta-LW green;
# grey = not significant),
# a "+"/"-" ionization suffix and the timepoint in each label, a solid vertical
# line at 0, and horizontal Q<0.05 (green) + P<0.05 (grey) threshold lines.
# (Previously this panel used coverage-% triangles with fill=ionization and
# colour=study x timepoint, which did not match the submission.)
# Data: our reproduced Table S5 (results/metaboanalyst/mummichog_s5/milk_mummichog_tableS5.csv).
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})
root <- here::here()
source(file.path(root, "figure-scripts/manuscript_figures/study_colors.R"))  # canonical study colours
source(file.path(root, "figure-scripts/0_figure-functions.R"))               # theme_imic() = Science-submission theme

d <- fread(file.path(root, "results/metaboanalyst/mummichog_s5/milk_mummichog_tableS5.csv"))
setnames(d, c("P-value", "Time Point", "Ionization Mode"), c("p_value", "timepoint", "ionization"))
d <- d[!is.na(`Enrichment Ratio`) & !is.na(p_value)]
d[, `:=`(
  fe       = `Enrichment Ratio`,                       # already signed by direction (neg = down)
  y        = -log10(p_value),
  study    = fifelse(Study == "MISAME", "MISAME-III", Study),
  ion_sign = fifelse(ionization == "Positive", "+", "-"),
  sig      = FDR < 0.05,
  # nominal (before-FDR) significance: raw P < 0.05 but not FDR-significant.
  nom_sig  = (p_value < 0.05) & !(FDR < 0.05)
)]
# Canonical study colours (ELICIT blue / MISAME orange / Mumta green); index by name.
studycols <- c(imic_study_cols, "Not significant" = "#CFCFCF")
# Colour ALL nominally-significant pathways (raw P < 0.05) by study, all filled
# (per Andrew, 2026-08-10) -- matching the submitted panel; the Q<0.05 (green) and
# P<0.05 (grey) reference lines mark which are FDR-significant. Non-significant -> grey.
d[, study_col := fcase(p_value < 0.05, study, default = "Not significant")]
d[, lab := fifelse(sig, paste0(Pathway, ion_sign, " (", timepoint, ")"), NA_character_)]
# Label the pathways DISCUSSED in the manuscript Results text (Fig 6B); each is
# FDR-significant here but was previously excluded by the old "top-3 by |fe|" rule.
# One label per distinct pathway (the most extreme |fe| instance). NB: "Prostaglandin
# formation from arachidonate" is FDR-significant only in the UP direction (fe > 0),
# contrary to the text's coordinated-reduction framing, so it is deliberately omitted.
discussed <- c("Arachidonic acid metabolism", "Linoleate metabolism",
               "De novo fatty acid biosynthesis", "Fatty acid activation",
               "Fatty Acid Metabolism", "Leukotriene metabolism")
lab_d <- d[sig == TRUE & Pathway %in% discussed][order(-abs(fe))][!duplicated(Pathway)]
q_line <- suppressWarnings(min(d[sig == TRUE, y]))
# Shared y-axis ceiling across Fig 6 Panels A/B/C, matching the submitted figure's style
# (all three panels share one -log10(P) height); see fig6A-untargeted-msea.R / fig6C-proteomics.R.
FIG6_Y_MAX <- 9

p <- ggplot(d, aes(fe, y, colour = study_col)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), colour = "#9A9A93") +
  geom_hline(yintercept = q_line, colour = "#3C8C3C") +
  annotate("text", x = min(d$fe), y = -log10(0.05), label = "P-value < 0.05",
           hjust = 0, vjust = -0.5, size = 2.1, colour = "#7A7A73") +
  annotate("text", x = min(d$fe), y = q_line, label = "Q-value < 0.05",
           hjust = 0, vjust = -0.5, size = 2.1, colour = "#3C8C3C") +
  geom_point(aes(shape = study_col), size = 2, alpha = 0.85) +   # study symbol = colourblind cue
  # size/label.padding trimmed 2026-08-26 (was 2.5/0.12) to better match Fig 3A/5A's
  # smaller, boxless repel-label text; base_size dropped 9->8 below for the same reason
  # (reported: Fig 6 labels read large next to the other figures).
  geom_label_repel(data = lab_d, aes(label = lab), size = 2.2, na.rm = TRUE,
                  max.overlaps = Inf, min.segment.length = 0, box.padding = 0.5,
                  label.padding = 0.1, label.size = 0.15, fill = "white",  # white box (matches Fig 6A)
                  force = 6, nudge_y = 0.6, segment.size = 0.25,
                  seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = studycols,
                      breaks = c("ELICIT", "MISAME-III", "Mumta-LW",
                                 "Not significant"),
                      name = "Study") +
  scale_shape_manual(values = imic_shapes_for(names(studycols)),
                     breaks = c("ELICIT", "MISAME-III", "Mumta-LW",
                                "Not significant"),
                     name = "Study") +
  scale_y_continuous(limits = c(0, FIG6_Y_MAX), breaks = seq(0, FIG6_Y_MAX, 1)) +
  labs(x = "Fold Enrichment", y = expression("–Log"[10]*"("*italic(P)*"-value)")) +
  theme_imic(base_size = 8) +   # Science-submission theme (Helvetica, font floors); matches Fig 3A/5A
  guides(colour = guide_legend(nrow = 2, byrow = TRUE),
         shape  = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.key.size = unit(0.35, "cm"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        plot.margin = margin(3, 6, 3, 3))

# Left-column half-page sub-panel of the full-page (7.25 in) Fig 6 (A/B/C left, D right).
ggsave(file.path(root, "figures/figure6_panelB_mummichog.png"),
       p, width = 210/2, height = 297/3, units = "mm", dpi = 300, bg = "white")  # 105x99 mm (submitted panel size)
cat("wrote figures/figure6_panelB_mummichog.png |", nrow(d), "pathways,",
    sum(d$sig), "FDR-sig,", nrow(lab_d), "labelled\n")
