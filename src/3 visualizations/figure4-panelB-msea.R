# figure4-panelB-msea.R
# =============================================================================
# SUPERSEDED 2026-09-08 -- do not run for the manuscript.
#   The live generator of Fig 3 Panel B is
#   figure-scripts/manuscript_figures/fig3B-pathway.R, which writes
#   figures/figure3_panelB_msea.png (renamed from the legacy
#   figure4_panelB_msea.png when the repo was aligned to the submitted
#   figure numbers). This script still writes the OLD name on purpose: it is
#   kept for reference only, and pointing it at the new name would give the
#   live panel a second, competing writer.
# =============================================================================
# =============================================================================
# Figure 4 / Panel B : "MSEA of B-vitamin pathways" (primary, combined arms).
#
# WHAT THIS IS
#   The metabolite-set-enrichment (MSEA / over-representation) panel of the
#   primary targeted B-vitamin / micronutrient pathways. It reproduces the
#   visualization STYLE of Trenton's manual figure
#   (`trenton scripts/old/imicPaperUntargetedMetabolomics.Rmd`, chunk
#   "*Untargeted Metabolites Figure"): a signed-enrichment-ratio volcano where
#     x = enrichment ratio (hits / expected), signed - for down-regulated cells,
#     y = -log10(raw p),
#   points are coloured by study, greyed when not nominally significant, the
#   nominal p = 0.05 line is drawn (red dashed) together with the enrichment-
#   ratio = 2 reference (blue dashed), and each nominally significant pathway is
#   labelled with its name.
#
# WHAT CHANGED (why this script exists)
#   Trenton built the panel by hand: he uploaded per-cell metabolite lists to
#   metaboanalyst.ca, downloaded the ORA result spreadsheets, and read them back
#   in the Rmd as `read_xlsx("../results/msea_ora_result_*.xlsx")` -- every one
#   of those files is flagged TODO-MISSING (Trenton-only intermediates never
#   committed to the repo, so the figure is not reproducible from this repo).
#   Our scripted MetaboAnalystR pipeline (`src/metaboanalyst/run-primary.R`) now
#   produces the identical ORA automatically, consolidated in
#   `results/metaboanalyst/primary_combined/primary_combined_supplementary_table.csv`.
#   This script renders the panel from THAT reproducible output instead of the
#   manual xlsx exports. It plots only what is in our CSVs -- no hand-entered
#   enrichment numbers.
#
# OUTPUT
#   figures/figure4_panelB_msea.png  -- sized ~1/2 A4 wide x 1/3 A4 tall
#   (105 x 99 mm), matching how panels B/C sit in the bottom row of
#   `src/3 visualizations/figure4-primary-volcano.R` compose_figure().
#
# Run from the repo root:
#   Rscript "src/3 visualizations/figure4-panelB-msea.R"
# =============================================================================

suppressMessages({
  library(dplyr); library(stringr); library(ggplot2); library(ggrepel); library(ggforce)
})

# ---------------------------------------------------------------------------
# Input: consolidated combined-arm primary ORA (one row per cell x pathway).
# Columns: study, studytime, direction, contrast, pathway, total, hits,
#          expected, raw_p, fdr_native, fdr_bh_pooled, significant, ...
# ---------------------------------------------------------------------------
supp <- read.csv(
  "results/metaboanalyst/primary_combined/primary_combined_supplementary_table.csv",
  check.names = FALSE, stringsAsFactors = FALSE)

# Tableau-20 palette (Trenton's), used to colour points by study.
tableau20 <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#86BCD6", "#FFBE7D", "#FF5850", "#A0CBE8", "#8CD17D",
  "#B6992D", "#499894", "#FABFD2", "#D37295", "#B7B7B7")

ALPHA <- 0.05          # nominal significance threshold
ER_REF <- 2            # enrichment-ratio reference (Trenton's blue dashed line)

# ---------------------------------------------------------------------------
# Prepare the plotting frame (mirrors Trenton's manual construction, but the
# input is our pipeline table rather than the msea_ora_result_*.xlsx files).
#   enrichment_ratio = hits / expected  (MetaboAnalyst's "Enrichment Ratio")
#   signed by direction so down-regulated cells sit on the negative x-axis.
# ---------------------------------------------------------------------------
set.seed(123)  # reproducible jitter for overlapping non-significant points

plot_df <- supp %>%
  mutate(
    study_label = recode(study,
                         Elicit = "ELICIT",
                         Misame = "MISAME-III",
                         Vital  = "Mumta-LW",
                         .default = study),
    enrichment_ratio  = hits / expected,
    enrichment_signed = if_else(tolower(direction) == "down",
                                -enrichment_ratio, enrichment_ratio),
    logP        = -log10(raw_p),
    is_sig      = raw_p < ALPHA,
    point_color = if_else(is_sig, study_label, "Not Significant"),
    # light jitter (Trenton's approach) so co-located non-significant pathways
    # -- many share the same hits/expected -- do not overplot into one dot.
    jitter_x = if_else(is_sig, enrichment_signed,
                       enrichment_signed + runif(n(), -0.6, 0.6)),
    jitter_y = if_else(is_sig, logP,
                       pmax(0, logP + runif(n(), -0.15, 0.15)))
  )

# Colour map: grey for non-significant, one tableau colour per study present.
studies_present <- sort(unique(plot_df$study_label[plot_df$is_sig]))
color_vals <- setNames(tableau20[seq_along(studies_present)], studies_present)
color_vals <- c("Not Significant" = "grey75", color_vals)

# Axis extents from the data (padded), keeping the ER = 2 reference visible.
x_hi <- max(c(plot_df$jitter_x, ER_REF + 1), na.rm = TRUE) * 1.05
x_lo <- min(c(plot_df$jitter_x, 0), na.rm = TRUE)
x_lo <- if (x_lo < 0) x_lo * 1.1 - 0.5 else -0.5
y_hi <- max(plot_df$logP, na.rm = TRUE) * 1.12

# Nominally significant pathways get a name label (Trenton labels hits>=1 &
# raw p < 0.05); collapse to one label per pathway x study to avoid duplicates
# from the same pathway recurring across visits.
label_df <- plot_df %>%
  filter(is_sig, hits >= 1) %>%
  group_by(study_label, pathway) %>%
  slice_max(logP, n = 1, with_ties = FALSE) %>%
  ungroup()

# Ellipses around "repeated" pathways (Trenton's MSEA/pathway figure feature):
# a pathway that is nominally significant (raw p < 0.05) in >= 2 analysis cells
# gets a black outline enclosing its significant points. One point per cell
# (study x visit x direction x contrast), so the count is the number of distinct
# significant appearances of that pathway.
repeated_paths <- plot_df %>%
  filter(is_sig) %>%
  count(pathway, name = "n_sig_cells") %>%
  filter(n_sig_cells >= 2) %>%
  pull(pathway)
ellipse_df <- plot_df %>% filter(is_sig, pathway %in% repeated_paths)
# geom_mark_ellipse needs >= 2 points per group; guard for the empty/all-unique case.
ellipse_layer <- if (nrow(ellipse_df) >= 2) {
  ggforce::geom_mark_ellipse(
    data = ellipse_df, aes(x = jitter_x, y = jitter_y, group = pathway),
    color = "black", fill = NA, linewidth = 0.35,
    expand = unit(2, "mm"), inherit.aes = FALSE, show.legend = FALSE)
} else NULL

# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------
panelB <- ggplot(plot_df, aes(x = jitter_x, y = jitter_y)) +
  ellipse_layer +
  geom_point(aes(color = point_color), size = 1.6, alpha = 0.85) +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", color = "red") +
  geom_vline(xintercept = ER_REF, linetype = "dashed", color = "blue") +
  geom_label_repel(
    data = label_df, aes(label = pathway, color = point_color),
    size = 2, label.padding = 0.12, box.padding = 0.4,
    min.segment.length = 0, max.overlaps = 200, show.legend = FALSE) +
  scale_color_manual(values = color_vals, name = "Study") +
  scale_x_continuous(limits = c(x_lo, x_hi)) +
  scale_y_continuous(limits = c(0, y_hi)) +
  labs(x = "Enrichment Ratio (signed by direction)",
       y = expression(-Log[10]*"(Raw P)")) +
  theme_classic(base_size = 7) +
  theme(
    axis.title   = element_text(size = 7),
    axis.text    = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    plot.margin  = margin(4, 6, 2, 4))

# ---------------------------------------------------------------------------
# Save: ~1/2 A4 width x 1/3 A4 height, UTF-8-safe ragg device, publication dpi.
# ---------------------------------------------------------------------------
ggsave(
  filename = "figures/figure4_panelB_msea.png",
  plot   = panelB,
  width  = 210 / 2,   # half A4 width  (mm)
  height = 297 / 3,   # one-third A4 height (mm)
  units  = "mm",
  dpi    = 600,
  device = ragg::agg_png)

cat("wrote figures/figure4_panelB_msea.png\n")
cat("pathways plotted:", nrow(plot_df),
    "| nominally significant (raw p <", ALPHA, "):", sum(plot_df$is_sig),
    "| FDR-significant:", sum(plot_df$significant %in% c(TRUE, "TRUE")), "\n")
