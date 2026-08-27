# figure-3b-pathway-submitted.R
# FAITHFUL port of Trenton's "Primary Outcomes (Pathway Analysis).Rmd" plot block
# (his 2026-08-03 script) -- the generator of the SUBMITTED Fig 3B. Ported verbatim:
# same study colours (#4E79A7 / #F28E2B / #E15759), theme_bw(7), the two annotated
# reference lines (P<0.05 grey, Q[FDR]<0.05 green), days->d / months->m timepoint
# abbreviation, ellipses (>=2 sig cells) under the points, and white boxed labels.
#
# Reads the same submitted MetaboAnalyst pathway exports already in the repo
# (results/metaboanalyst/primary_pathway_trenton/*.csv) that fig3B-pathway.R
# uses, so the DATA is identical -- only the styling is Trenton's exact original.
# Output: figures/trenton_ports/figure_3b_pathway_submitted.png (+ comparison crop).
suppressMessages({ library(dplyr); library(ggplot2); library(ggrepel); library(ggforce); library(stringr); library(readr) })

DIR <- "results/metaboanalyst/primary_pathway_trenton"
OUT1 <- "figures/trenton_ports/figure_3b_pathway_submitted.png"
OUT2 <- "Manuscript/figure_comparison/msea_compare/fig3B_faithful.png"

.parse_cell <- function(fname) {
  study_key <- str_extract(fname, "^(Elicit|Misame|Mumpta)")
  study <- dplyr::recode(study_key, Elicit = "ELICIT", Misame = "MISAME-III", Mumpta = "Mumta-LW")
  tp <- fname %>% str_remove("^(Elicit|Misame|Mumpta)_") %>% str_remove("_Combined_?\\.csv$")
  tp <- dplyr::recode(tp,
    "1_month" = "1 month", "5_months" = "5 months", "14_21_days" = "14-21 days",
    "1_2_months" = "1-2 months", "3_4_months" = "3-4 months",
    "1_5_months" = "1.5 months", "2_months" = "2 months", .default = tp)
  list(study = study, timepoint = tp)
}

files <- list.files(DIR, pattern = "\\.csv$", full.names = TRUE)
primary_pathanaly_results <- bind_rows(lapply(files, function(f) {
  x <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE); names(x)[1] <- "pathway"
  ce <- .parse_cell(basename(f))
  tibble(study = ce$study, timepoint = ce$timepoint, pathway = x$pathway,
         raw_p = as.numeric(x[["Raw p"]]), fdr = as.numeric(x[["FDR"]]),
         impact = as.numeric(x[["Impact"]]), logp = -log10(as.numeric(x[["Raw p"]])))
}))

# ---- Trenton's exact styling block (Primary Outcomes (Pathway Analysis).Rmd) ----
study_colors <- c("ELICIT" = "#4E79A7", "MISAME-III" = "#F28E2B", "Mumta-LW" = "#E15759")
p_line_col <- "#BAB0AC"; q_line_col <- "#59A14F"
p_thr <- 0.05; y_p_line <- -log10(p_thr)
fdr_p_thr <- primary_pathanaly_results %>% filter(fdr < 0.05) %>%
  summarise(v = max(raw_p, na.rm = TRUE)) %>% pull(v)
y_fdr_line <- -log10(fdr_p_thr)

plot_df <- primary_pathanaly_results %>%
  filter(!is.na(raw_p), !is.na(impact)) %>%
  mutate(
    study = factor(study, levels = c("ELICIT", "MISAME-III", "Mumta-LW")),
    timepoint_clean = timepoint %>% str_replace_all("days", "d") %>%
      str_replace_all("months", "m") %>% str_replace_all("month", "m") %>% str_squish(),
    pathway_short = recode(pathway,
      "Nicotinate and nicotinamide metabolism" = "NAD/NAM Met",
      "Riboflavin metabolism" = "Riboflavin Met",
      "Thiamine metabolism" = "Thiamine Met",
      "Vitamin B6 metabolism" = "Vit B6 Met", .default = pathway),
    pathway_label = if_else(raw_p < p_thr, paste0(pathway_short, " (", timepoint_clean, ")"), NA_character_),
    color_group = if_else(raw_p < p_thr, as.character(study), "Not significant"))

ellipse_df <- plot_df %>% filter(raw_p < p_thr) %>%
  add_count(pathway, name = "pc") %>% filter(pc >= 2)
full_colors <- c("Not significant" = p_line_col, study_colors)

# RECONSTRUCT the submitted 3B look (the exact May pathway script is missing; his
# Aug Primary Outcomes (Pathway Analysis).Rmd drifted). Three deltas restored:
#  (1) integer y-axis ticks; (2) green line labelled "P-value < <raw-p thr>" (his
#  Aug relabel was "Q[FDR] < 0.05"); (3) font scaled by s = 105/240 so his 240 mm
#  base_size=7 proportions render correctly at the 105 mm submission panel size.
s <- 105 / 240
p <- ggplot(plot_df, aes(x = impact, y = logp, color = color_group)) +
  geom_hline(yintercept = y_p_line, color = p_line_col, linewidth = 0.5 * s) +
  geom_hline(yintercept = y_fdr_line, color = q_line_col, linewidth = 0.5 * s) +
  annotate("text", x = 0, y = y_p_line + 0.12, label = "italic(P)*\"-value\" < 0.05",
           parse = TRUE, hjust = 0, size = 2 * s, color = p_line_col) +
  annotate("text", x = 0, y = y_fdr_line + 0.12,
           label = paste0("italic(P)*\"-value\" < ", signif(fdr_p_thr, 2)),
           parse = TRUE, hjust = 0, size = 2 * s, color = q_line_col) +
  ggforce::geom_mark_ellipse(data = ellipse_df, aes(group = pathway), color = "black",
                             fill = NA, linewidth = 0.35 * s, expand = grid::unit(2, "mm"),
                             inherit.aes = TRUE, show.legend = FALSE) +
  geom_point(size = 1.8 * s, alpha = 0.9) +
  geom_label_repel(data = ellipse_df, aes(label = pathway_label, color = study),
                   size = 2 * s, fill = "white", label.size = 0.2 * s, box.padding = 0.35,
                   point.padding = 0.3, segment.size = 0.4 * s, min.segment.length = 0,
                   max.overlaps = Inf, seed = 123, show.legend = FALSE) +
  scale_color_manual(values = full_colors) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.05, 1.05)) +
  scale_y_continuous(breaks = seq(1, 11, 1)) +
  labs(x = "Pathway Impact", y = expression(-log[10](italic(P)*"-value")), color = "Study") +
  theme_bw(base_size = 7 * s) +
  theme(axis.text = element_text(size = 7 * s), axis.title = element_text(size = 7 * s),
        legend.text = element_text(size = 7 * s), legend.title = element_text(size = 7 * s),
        legend.position = "bottom", legend.box = "horizontal", plot.margin = margin(5, 5, 5, 5)) +
  guides(color = guide_legend(override.aes = list(size = 3 * s)))

dir.create(dirname(OUT1), recursive = TRUE, showWarnings = FALSE)
# Panel dimensions from Trenton's SUBMISSION-ERA panel convention, confirmed in code
# in imicPaperTriglycerides.Rmd and imicPaperProteomics.Rmd (both dated the May 26
# submission): ggsave(width = 210/2, height = 297/3, units = "mm", dpi = 600) --
# half-A4 width x one-third-A4 height (105 x 99 mm, aspect 1.06, near-square).
# NOTE: no pathway-specific export script exists (the only pathway-plot code, the
# Aug "Primary Outcomes (Pathway Analysis).Rmd", uses 240x135 mm = 1.778, the wide
# version); we adopt his documented panel convention so 3B is sized like the others.
for (o in c(OUT1, OUT2))
  ggsave(o, p, width = 210/2, height = 297/3, units = "mm", dpi = 600, bg = "white")
cat("wrote faithful 3B (Trenton panel convention 105x99mm) | cells:", length(files),
    "| ellipse groups:", dplyr::n_distinct(ellipse_df$pathway),
    "| green line raw-p:", signif(fdr_p_thr, 3), "\n")
