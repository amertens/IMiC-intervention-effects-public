# =============================================================================
# figS3-hmo-bioactives.R
#
# Builds Fig S3: forest plots of the combined-arm intervention effects (Z-scored
# outcomes, 95% CI) on the secondary outcomes, one facet per study and visit.
# Panel A: human milk oligosaccharides; Panel B: bioactive proteins. Points are
# coloured by significance as in Fig 2; no secondary outcome is FDR-significant.
#
# Inputs:  results/subsetted results/secondary_hmo.csv
#          results/subsetted results/secondary_bioactives.csv
# Outputs: figures/figureS3_hmo_bioactives.{png,pdf,eps}
# =============================================================================

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
library(cowplot)

# Standalone significance legend, built from a tiny synthetic data frame rather than
# taken from either forest plot's own legend: in ggplot2 4.0.x a discrete level with
# zero data rows (here "Sig", since nothing in this figure clears FDR) draws no key
# glyph even with drop = FALSE. Same pattern as create_study_legend() in
# fig5-tertiary-composite.R. Key and facet titles are worded as in Fig 2.
SIG_LABELS <- c("Not Significant" = "Not significant", "Sig before FDR" = "Nominally significant",
                "Sig" = "FDR-significant")

create_significance_legend <- function() {
  lev <- c("Not Significant", "Sig before FDR", "Sig")
  legend_df <- data.frame(x = 1, y = 1:3, sigcat = factor(lev, levels = lev))
  p <- ggplot(legend_df, aes(x, y, colour = sigcat, shape = sigcat)) +
    geom_point(size = 2) +
    scale_shape_manual(values = c("Not Significant" = 1, "Sig before FDR" = 1, "Sig" = 19),
                       labels = SIG_LABELS, name = "Statistical Significance") +
    scale_colour_manual(values = c("Not Significant" = "grey60", "Sig before FDR" = tableau10[1],
                                   "Sig" = tableau10[2]),
                        labels = SIG_LABELS, name = "Statistical Significance") +
    theme_void(base_size = 8) +
    theme(legend.position = "right", legend.title = element_text(size = 8),
          legend.text = element_text(size = 7))
  comps <- cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
  good  <- Filter(function(g) !inherits(g, "zeroGrob"), comps)
  if (length(good) == 0) stop("create_significance_legend(): no non-empty guide-box found")
  good[[1]]
}

forest_plot <- function(tab, pos="none", y_limits=NULL){

  # canonical biomarker labels (0_figure-functions.R)
  tab$label_f <- canonical_label(tab$biomarker, tab$label_f)
  tab$sigFDR <- factor(tab$sigFDR, levels=c(0,1))
  tab$sig <- factor(tab$sig, levels=c(0,1))
  tab$study <- dplyr::recode(as.character(tab$study),
                             "Vital"="Mumta-LW", "Elicit"="ELICIT", "Misame"="MISAME-III")
  tab$study <- factor(tab$study, levels=c("Mumta-LW","ELICIT","MISAME-III"))
  tab$visit <- factor(tab$visit, levels=c(  "1.5 mo.","2 mo.","1 mo.","5 mo.","14-21 days", "1-2 mo.","3-4 mo."))
  # facet titles written "Study (Timepoint)", as in Fig 2
  tab <- tab %>% arrange(study, visit) %>%
    mutate(panel = factor(paste0(study, " (", visit, ")"),
                          levels = unique(paste0(study, " (", visit, ")"))))
  tab=tab %>% mutate(sigcat=case_when(
    sigFDR == 1 & sig == 1 ~ "Sig",
    sigFDR == 0 & sig == 1 ~ "Sig before FDR",
    sigFDR == 0 & sig == 0 ~ "Not Significant"
  ), sigcat=factor(sigcat, levels=c("Not Significant", "Sig before FDR", "Sig")))

  # 3-tier scheme consistent with Fig 2 and Fig S2: grey = not significant,
  # tableau blue = significant before FDR, tableau orange = significant after FDR
  # (open circle before FDR, filled after). theme_bw() draws a border around each
  # facet.
  ggplot(tab, aes(x = reorder(label_f, -est), y = est, color=sigcat, shape=sigcat, group=contrast)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = cil, ymax = ciu), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    coord_flip() +
    scale_shape_manual(values = c("Not Significant"=1, "Sig before FDR"=1, "Sig"=19),
                       name="Statistical Significance", labels = SIG_LABELS, drop = FALSE,
                       limits = c("Not Significant","Sig before FDR","Sig")) +  # force orange "Sig" key even with no sig data
    scale_color_manual(values = c("Not Significant"="grey60",
                                  "Sig before FDR"=tableau10[1], "Sig"=tableau10[2]),
                       name="Statistical Significance", labels = SIG_LABELS, drop = FALSE,
                       limits = c("Not Significant","Sig before FDR","Sig")) +
    facet_wrap(~ panel, ncol = 4) +
    { if (!is.null(y_limits)) scale_y_continuous(limits = y_limits) } +  # shared x-axis (post-flip) across A/B
    labs(x = "", y = "") +
    # no horizontal gridlines or y ticks, as Fig 2
    theme_bw() + theme(legend.position = pos, axis.text.y=element_text(size=7),
                       panel.grid.major.y = element_blank(), axis.ticks.y = element_blank(),
                       panel.grid.minor = element_blank(),
                       strip.background = element_blank(),   # white strips, as Fig 2
                       strip.text = element_text(size = 7),
                       legend.title = element_text(size = 8),
                       legend.text = element_text(size = 7))
}


hmo_dat <- read.csv(here("results/subsetted results/secondary_hmo.csv")) %>%
  dplyr::filter(biomarker != "secretor")
bioact_dat <- read.csv(here("results/subsetted results/secondary_bioactives.csv"))
# Shared x-axis (post-coord_flip) range across Panels A and B, so both read on the
# same scale instead of each auto-scaling to its own data. No manual padding:
# scale_y_continuous()'s default expansion adds it.
shared_ylim <- range(c(hmo_dat$cil, hmo_dat$ciu, bioact_dat$cil, bioact_dat$ciu), na.rm = TRUE)

# Panel A: HMOs. "Maternal Secretor Status" (biomarker "secretor") is dropped above,
# as in the submitted figure.
# pos="none": both panels' own legends are suppressed; the figure shows one legend,
# the hand-built create_significance_legend() composited in below.
p4= forest_plot(hmo_dat, pos="none", y_limits = shared_ylim)

# Panel B: bioactive proteins
p5= forest_plot(bioact_dat, pos="none", y_limits = shared_ylim)

# ---- Fig S3: secondary outcomes (HMOs + bioactive proteins), combined arms ----
# p4 and p5 are aligned as two raw ggplots (align="v", axis="lr"); cowplot's align
# step needs raw ggplot cells, so the legend is not nested into this grid.
stack_ab <- plot_grid(
  p4+ theme(plot.margin = margin(0,0,0,0)),
  p5+ theme(plot.margin = margin(0,0,0,0)),
  ncol = 1,
  labels = "AUTO",
  rel_heights = c(1,0.4),
  align = "v",           # Align vertically
  axis = "lr"            # Align left and right axes
)
# The legend is overlaid after alignment, on the empty facet cell in Panel A's 4-column
# grid (row 2 has only MISAME-III's three visits, leaving column 4 empty).
fig_s3 <- ggdraw(stack_ab) +
  draw_plot(create_significance_legend(), x = 0.78, y = 0.30, width = 0.20, height = 0.22)

save_figure_3way(fig_s3, name = "figureS3_hmo_bioactives", width = 7.25, height = 10,
                 dir = paste0(here::here(), "/figures"))
