# =============================================================================
# fig2-primary-forest.R
#
# Builds Fig 2: forest plots of the average treatment effects (combined
# intervention arms vs control, Z-scored outcomes, 95% CI) on the primary milk
# outcomes, one facet per study and collection visit. Panel A: macronutrients;
# B: micronutrients and fat-soluble vitamins; C: B-vitamins and their metabolites.
# Points are coloured by significance (grey = not significant, blue = P < 0.05
# before FDR, orange = FDR-significant).
#
# Inputs:  results/subsetted results/primary_{macro,micro,bvit}.csv
# Outputs: figures/figure2.{png,pdf,eps}
# =============================================================================

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
library(cowplot)

SIG_LABELS <- c("Not Significant" = "Not significant", "Sig before FDR" = "Nominally significant",
                "Sig" = "FDR-significant")

forest_plot <- function(tab){

  # canonical biomarker labels (0_figure-functions.R)
  tab$label_f <- canonical_label(tab$biomarker, tab$label_f)
  tab$sigFDR <- factor(tab$sigFDR, levels=c(0,1))
  tab$sig <- factor(tab$sig, levels=c(0,1))
  # published study names (the raw keys are Vital/Elicit/Misame)
  tab$study <- dplyr::recode(as.character(tab$study),
                             "Vital"="Mumta-LW", "Elicit"="ELICIT", "Misame"="MISAME-III")
  tab$study <- factor(tab$study, levels=c("Mumta-LW","ELICIT","MISAME-III"))
  tab$visit <- factor(tab$visit, levels=c(  "1.5 mo.","2 mo.","1 mo.","5 mo.","14-21 days", "1-2 mo.","3-4 mo."))
  # facet titles written "Study (Timepoint)", as in Fig 1
  tab <- tab %>% arrange(study, visit) %>%
    mutate(panel = factor(paste0(study, " (", visit, ")"),
                          levels = unique(paste0(study, " (", visit, ")"))))
  tab=tab %>% mutate(sigcat=case_when(
    sigFDR == 1 & sig == 1 ~ "Sig",
    sigFDR == 0 & sig == 1 ~ "Sig before FDR",
    sigFDR == 0 & sig == 0 ~ "Not Significant"
  ), sigcat=factor(sigcat, levels=c("Not Significant", "Sig before FDR", "Sig")))

  # 3-tier colouring (matches the paper's volcanoes): grey = not significant,
  # blue = significant before FDR, orange = significant after FDR.
  ggplot(tab, aes(x = reorder(label_f, -est), y = est, color=sigcat, shape=sigcat, group=contrast)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = cil, ymax = ciu), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    coord_flip() +
    scale_shape_manual(values = c("Not Significant"=1, "Sig before FDR"=1, "Sig"=19),
                       labels = SIG_LABELS, drop = FALSE) +   # open circle for before-FDR
    scale_color_manual(values = c("Not Significant"="grey60",
                                  "Sig before FDR"=tableau10[1],
                                  "Sig"=tableau10[2]), labels = SIG_LABELS, drop = FALSE) +
    scale_x_discrete(labels = plotmath_expr) +   # B-vitamin subscripts via plotmath (Arial has no subscript-digit glyphs for PDF/EPS)
    facet_wrap(~ panel, ncol = 4) +
    labs(x = "",
         y = "") +
    # no horizontal gridlines or y ticks, as in the other figures
    theme_bw() + theme(legend.position = "none", axis.text.y=element_text(size=7),
                       panel.grid.major.y = element_blank(), axis.ticks.y = element_blank(),
                       panel.grid.minor = element_blank(),
                       strip.background = element_blank(),      # white strips
                       strip.text = element_text(size = 7))     # theme_bw() = border around each facet
}


# Panel A: macronutrients
p1= forest_plot(read.csv(here("results/subsetted results/primary_macro.csv")))

# Panel B: micronutrients and fat-soluble vitamins
p2= forest_plot(read.csv(here("results/subsetted results/primary_micro.csv")))

# Panel C: B-vitamins
p3= forest_plot(read.csv(here("results/subsetted results/primary_bvit.csv")))


# Shared significance legend placed inside the B-vitamin panel's empty bottom-
# right facet cell (the MISAME-III row has only 3 of 4 columns) rather than in a
# separate strip.
p3_leg <- p3 +
  theme(legend.position = c(0.995, 0.18),        # low in the empty bottom-right cell, clear of the last panel's plotted points
        legend.justification = c(1, 0.5),        # right-anchored so it can't spill left
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank(),
        legend.key.size = unit(0.35, "cm"),       # keys sized to match the panel's own axis/strip text, not larger
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  guides(colour = guide_legend(title = "Statistical Significance"),
         shape  = guide_legend(title = "Statistical Significance"))

fig2 <- plot_grid(
  p1 + theme(plot.margin = margin(0,6,-1,0)),
  p2 + theme(plot.margin = margin(0,6,-1,0)),
  p3_leg + theme(plot.margin = margin(0,6,0,0)),
  ncol = 1,
  labels = "AUTO",
  # Panel C (22 B-vitamin rows) gets the most height so its row labels don't
  # collide.
  rel_heights = c(0.30, 0.66, 0.95),
  align = "v",           # Align vertically
  axis = "lr")           # Align left and right axes

# Export PDF + EPS + PNG (white background).
save_figure_3way(fig2, name = "figure2", width = 7.25, height = 10.6)   # taller for Panel C
