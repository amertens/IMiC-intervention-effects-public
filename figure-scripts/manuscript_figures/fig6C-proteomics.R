# fig6C-proteomics.R
# Fig 6 Panel C: milk-proteome GO Biological-Process over-representation, VOLCANO.
# Style matches the SUBMITTED Fig 6C and Panels A/B (x = signed fold enrichment,
# y = -log10(P), colour = study, P<0.05 + Q<0.05 threshold lines, extremes labelled).
# Data = UniProt-native rerun (results/proteomics_go_uniprot.csv, from
# src/2 analysis/55-proteomics-go-uniprot.R), Trenton's finished method: each protein
# counted once. NOTE vs submitted: under this stricter method the up-regulated
# anti-proteolysis / immune terms no longer clear FDR; only down-regulated terms
# survive (MISAME telomere/DNA, Mumta nucleotide/NAD/energy).
# Out: figures/figure6_panelC_proteomics_go.png
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})
root <- here::here()
source(file.path(root, "figure-scripts/manuscript_figures/study_colors.R"))  # canonical study colours
source(file.path(root, "figure-scripts/0_figure-functions.R"))               # theme_imic() = Science-submission theme

d <- fread(file.path(root, "results/proteomics_go_uniprot.csv"))
d[, `:=`(fe = signed_fold_enrichment, y = -log10(pvalue), sig = p.adjust < 0.05)]
# study colours: canonical map (match graphical abstract). Only MISAME-III + Mumta-LW
# appear here; index BY NAME so each keeps its colour (orange / purple), never by order.
studycols <- imic_study_cols[c("MISAME-III", "Mumta-LW")]  # canonical map (ELICIT absent from proteomics)
# SAME THREE TIERS AS PANEL A (2026-09-25, per Andrew), so the one shared key under
# Panel C is exactly right for A-C: FDR-significant (Q < 0.05) -> study colour + study
# symbol; significant before FDR only (P < 0.05) -> open light-blue circle; otherwise grey.
# (Until then every P < 0.05 point was filled in study colour, per Andrew 2026-08-10, to
# match the submitted panel -- which the shared Panel A key would have misdescribed.)
d[, study_col := fcase(sig, study, pvalue < 0.05, "Sig before FDR", default = "Not Significant")]
setorder(d[, draw := fcase(study_col == "Not Significant", 0L,
                           study_col == "Sig before FDR", 1L, default = 2L)], draw)  # filled on top
q_line <- suppressWarnings(min(d[sig == TRUE, y]))     # -log10(P) at the BH-FDR 0.05 frontier
# Trenton's abbr_proteomics() (from Exploratory Outcomes (Proteomics - UniProt).Rmd),
# ported verbatim so the GO labels read as his current panel.
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
# Label FDR-significant terms: top 3 per study x direction by |fold enrichment| (fewer
# than his standalone's 10/direction, since this is the half-column composite panel).
labs <- d[sig == TRUE][order(-abs(fe))][, .SD[!duplicated(lab)], by = .(study, direction)][
  , head(.SD, 3), by = .(study, direction)]

# Shared y-axis ceiling across Fig 6 Panels A/B/C, matching the submitted figure's style
# (all three panels share one -log10(P) height); see fig6A-untargeted-msea.R / fig6B-mummichog.R.
FIG6_Y_MAX <- 9
# No legend of its own: one shared Fig 6 A-C key (Panel A's) sits under this panel, and
# A/B/C share one panel size; see fig6A-untargeted-msea.R (2026-09-25).
FIG6_PANEL_H_MM <- 95
pal <- c(studycols, "Sig before FDR" = "#6BAED6", "Not Significant" = "grey75")
tier_shapes <- imic_shapes_for(names(pal))
tier_shapes["Sig before FDR"] <- 1   # open circle, as in Panel A
p <- ggplot(d, aes(fe, y)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "solid", colour = "#9A9A93") +
  geom_hline(yintercept = q_line, linetype = "solid", colour = "#3C8C3C") +
  annotate("text", x = min(d$fe), y = -log10(0.05), label = "P-value < 0.05",
           hjust = 0, vjust = -0.5, size = 2.1, colour = "#7A7A73") +
  annotate("text", x = min(d$fe), y = q_line, label = "Q-value < 0.05",
           hjust = 0, vjust = -0.5, size = 2.1, colour = "#3C8C3C") +
  geom_point(aes(colour = study_col, shape = study_col), size = 2, alpha = 0.85) +  # study symbol = colourblind cue
  # size/label.padding trimmed 2026-08-26 (was 2.5/0.12), matching fig6B-mummichog.R's
  # rationale -- Fig 3A/5A's smaller, boxless repel labels; base_size dropped 9->8 below.
  geom_label_repel(data = labs, aes(label = lab, colour = study), size = 2.2,
                  max.overlaps = Inf, min.segment.length = 0, box.padding = 0.3,
                  label.padding = 0.1, label.size = 0.15, fill = "white",  # white box (matches Fig 6A)
                  seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = pal,
                      breaks = c("MISAME-III", "Mumta-LW", "Sig before FDR", "Not Significant"),
                      name = "Study") +
  scale_shape_manual(values = tier_shapes,
                     breaks = c("MISAME-III", "Mumta-LW", "Sig before FDR", "Not Significant"),
                     name = "Study") +
  scale_y_continuous(limits = c(0, FIG6_Y_MAX), breaks = seq(0, FIG6_Y_MAX, 1)) +  # shared Fig 6 A/B/C height
  labs(x = "Fold Enrichment (signed by direction)",
       y = expression("–Log"[10]*"("*italic(P)*"-value)")) +
  theme_imic(base_size = 8) +   # Science-submission theme (Helvetica, font floors); matches Fig 3A/5A
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "none", legend.key.size = unit(0.35, "cm"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        plot.margin = margin(4, 6, 4, 4))

# Left-column half-page sub-panel of the full-page (7.25 in) Fig 6 (A/B/C left, D right).
ggsave(file.path(root, "figures/figure6_panelC_proteomics_go.png"),
       p, width = 210/2, height = FIG6_PANEL_H_MM, units = "mm", dpi = 300, bg = "white",  # 105 mm wide, shared A/B/C height
       device = ragg::agg_png)
cat("wrote figures/figure6_panelC_proteomics_go.png |", nrow(d), "terms,",
    sum(d$sig), "FDR-sig,", nrow(labs), "labelled | up FDR-sig:",
    d[sig == TRUE & direction == "Upregulated", .N], "\n")
