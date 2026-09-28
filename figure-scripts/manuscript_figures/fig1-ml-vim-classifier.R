# =============================================================================
# fig1-ml-vim-classifier.R
#
# Builds Fig 1. Panel A: cross-validated AUC (95% CI) of Super Learner classifiers
# that predict trial arm from each group of milk analytes, by study, contrast and
# collection time; ELICIT includes the infant-azithromycin arm as a negative control.
# Panel B: adjusted intervention effect (ATE, 95% CI) on the first principal
# component of each milk modality.
#
# Inputs:  figure-data/SL_vim_plot_data.RDS (written by src/3 visualizations/5-SL_VIM_plots.R)
#          results/pca_intervention_effects_results.RDS (src/2 analysis/3_adjusted_analysis_pca.R)
# Outputs: figures/figure1.{png,pdf,eps}
# [needs restricted data] SL_vim_plot_data.RDS derives from
# results/SL_individual_lab_vim_res.RDS, which src/2 analysis/3b-SL_vim_individual_lab.R
# fits on participant-level data and which is not shipped.
# =============================================================================
suppressPackageStartupMessages(library(here))
# the relative '../../' paths below resolve from this folder
setwd(file.path(here::here(), "figure-scripts/manuscript_figures"))

library(tidyverse)
library(cowplot)

# shared palette (tableau10), theme_imic(), science_dims and save_figure_3way()
source(file.path(here::here(), "figure-scripts/0_figure-functions.R"))

# Science 2-column full-page width
page_w <- science_dims$full_page$width        # 7.25 in (184 mm)

# SL variable-importance plot data, including the ELICIT infant-azithromycin
# negative-control arm.
SLvim_all <- readRDS("../../figure-data/SL_vim_plot_data.RDS")
# "Pre+Postnatal" -> "Pre + Postnatal", matching the spaced "+" of the other
# combined arms
levels(SLvim_all$arm_f) <- sub("Pre+Postnatal", "Pre + Postnatal", levels(SLvim_all$arm_f), fixed = TRUE)

# Shared CV-AUC x-range that covers EVERY confidence interval across the three
# studies (some ci.lb reach ~0.26), so nothing is clipped. Data-driven with
# minimal padding to keep white space down.
slvim_xlim <- c(floor(min(SLvim_all$ci.lb, na.rm = TRUE) / 0.05) * 0.05,
                min(1.0, ceiling(max(SLvim_all$ci.ub, na.rm = TRUE) / 0.05) * 0.05))

# ---- Panel A -----------------------------------------------------------------
# Rebuilt from the plot data rather than from the stored ggplot objects, which
# were made with an older ggplot2 and fail grid.draw under the current version.
build_slvim_panel <- function(d, title = "", show_x = FALSE,
                              legend_pos = "none", xlab_text = "",
                              xlim_range = slvim_xlim) {
  # Order predictor groups so "All" stays at the top of the y-axis
  group_levels <- c("All", "Macronutrients", "Micronutrients",
                    "B-vitamins", "HMOs", "Targeted proteins",
                    "Targeted metabolomics")
  d$group <- factor(d$group, levels = rev(intersect(group_levels, unique(d$group))))
  ggplot(d, aes(x = cvAUC, y = group, colour = visit_f)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "black", linewidth = 0.5) +  # same as Panel B
    geom_linerange(aes(xmin = ci.lb, xmax = ci.ub),
                   position = position_dodge(width = 0.6)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
    facet_grid(. ~ arm_f) +
    scale_colour_manual(
      values = c("<1 month" = "#7F7F7F",     # grey / orange / teal
                 "1-2 months" = "#FF7F0E",
                 "2-5 months" = "#17BECF"),
      name = "Collection time") +
    scale_x_continuous(breaks = c(0.4, 0.6, 0.8, 1.0),
                       expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(xlim = xlim_range) +   # zoom, not filter -> keeps every CI
    labs(x = if (show_x) "CV-AUC" else xlab_text, y = title) +
    # white strips (Fig 2's forest style); no horizontal/minor gridlines or y ticks,
    # as in the other figures
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          legend.position  = legend_pos,
          legend.text      = element_text(size = 7),
          legend.title     = element_text(size = 7))
}

# ---- Panel B -----------------------------------------------------------------
# ATEs on the first principal component of all 8 milk modalities (including
# untargeted metabolomics and microbiome).
pca_df <- readRDS(file="../../results/pca_intervention_effects_results.RDS") %>% filter(measure=="ATE")
pca_df$studytime <- gsub("-40", "-1.5", pca_df$studytime)
pca_df$studytime <- gsub("-56", "-2", pca_df$studytime)
pca_df$studytime <- gsub("Vital", "Mumta-LW", pca_df$studytime)
pca_df$studytime <- gsub("Misame-1", "Misame (14-21 days)", pca_df$studytime)
pca_df$studytime <- gsub("Misame-2", "Misame (1-2 mo.)", pca_df$studytime)
pca_df$studytime <- gsub("Misame-3", "Misame (3-4 mo.)", pca_df$studytime)
pca_df$studytime <- gsub("-1.5", " (1.5 mo.)", pca_df$studytime)
pca_df$studytime <- gsub("Elicit-1", "ELICIT (1 mo.)", pca_df$studytime)
pca_df$studytime <- gsub("Mumta-LW-2", "Mumta-LW (2 mo.)", pca_df$studytime)
pca_df$studytime <- gsub("Elicit-5", "ELICIT (5 mo.)", pca_df$studytime)

p_elicit <- build_slvim_panel(dplyr::filter(SLvim_all, studyid == "ELICIT"),    title = "ELICIT")
p_vital  <- build_slvim_panel(dplyr::filter(SLvim_all, studyid == "Mumta-LW"),  title = "Mumta-LW")
p_misame <- build_slvim_panel(dplyr::filter(SLvim_all, studyid == "Misame-III"), title = "MISAME-III",
                              show_x = TRUE, legend_pos = "bottom")

SLvim_lab_plot <- plot_grid(p_elicit, p_vital, p_misame,
                            ncol = 1, rel_heights = c(1, 1, 1.25),
                            align = "v", axis = "lr")

pca_df <- pca_df %>% mutate(
  studytime=gsub("Misame","MISAME-III",studytime),
  studyid=case_when(
    grepl("ELICIT",studytime) ~"ELICIT",
    grepl("Mumta-LW",studytime) ~"Mumta-LW",
    grepl("MISAME",studytime) ~"MISAME-III"
  )
  )


# Map raw biomarker keys ("bvit_pca", "macro_pca", ...) to user-facing
# category labels expected by the factor levels below. Without this remap,
# factor() set every row to NA and Panel B lost its category axis + coloring.
.lab_map <- c(macro_pca="Macronutrients", micro_pca="Micronutrients", bvit_pca="B-vitamins",
              HMO_pca="HMOs", protein_pca="Proteins", metabolomics_pca="Targeted metabolomics",
              untarget_metabolomics_pca="Untargeted metabolomics", untargeted_pca="Untargeted metabolomics",
              microbiome_pca="Microbiome")
# Map raw keys -> user-facing labels; values already in final form (from the
# results/ copy) are not in the map and are kept via coalesce.
pca_df <- pca_df %>% mutate(
  label_f = dplyr::coalesce(unname(.lab_map[as.character(label_f)]), as.character(label_f)))

#add blank rows to vital and elicit for blank facet
df_blank <- data.frame(studyid=c("ELICIT","Mumta-LW"),
  studytime=c("",""), label_f=c("Macronutrients","Macronutrients"))
pca_df <- bind_rows(pca_df, df_blank)
pca_df <- pca_df %>% mutate(
  label_f=factor(label_f, levels=rev(c( "Macronutrients","Micronutrients","B-vitamins", "HMOs","Proteins", "Targeted metabolomics" ,"Untargeted metabolomics","Microbiome"))))


p_pca1 <- ggplot(pca_df %>% filter(studyid=="ELICIT"), aes(x=label_f, y=est, color=label_f)) + geom_point() +
  geom_linerange(aes(ymin=cil, ymax=ciu)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  coord_flip() +
  facet_grid(~studytime) +
  scale_color_manual(values=rev(tableau10)) +
  theme_bw(base_size = 8, base_family = "Helvetica") +   # same theme as Panel A
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.position = "none") +
  ylab("") + xlab("")


p_pca2 <- ggplot(pca_df %>% filter(studyid=="Mumta-LW"), aes(x=label_f, y=est, color=label_f)) + geom_point() +
  geom_linerange(aes(ymin=cil, ymax=ciu)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  coord_flip() +
  facet_grid(~studytime) +
  scale_color_manual(values=rev(tableau10)) +
  theme_bw(base_size = 8, base_family = "Helvetica") +   # same theme as Panel A
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.position = "none") +
  ylab("") + xlab("Milk Modality")

p_pca3 <- ggplot(pca_df %>% filter(studyid=="MISAME-III") %>%
                   mutate(studytime=factor(studytime, levels = c("MISAME-III (14-21 days)","MISAME-III (1-2 mo.)","MISAME-III (3-4 mo.)" ))), aes(x=label_f, y=est, color=label_f)) + geom_point() +
  geom_linerange(aes(ymin=cil, ymax=ciu)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  coord_flip() +
  facet_grid(~studytime) +
  scale_color_manual(values=rev(tableau10)) +
  theme_bw(base_size = 8, base_family = "Helvetica") +   # same theme as Panel A
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.position = "none") +
  ylab("Average Treatment Effect") + xlab("")


p_pca <- plot_grid(p_pca1 ,
                            p_pca2 , 
                            p_pca3,
                            ncol=1, rel_heights=c(1,1,1.1))

fig_1 <- plot_grid(
  SLvim_lab_plot,
  p_pca,
  align="h",
  axis ="l",
      labels="AUTO",
  ncol = 1, nrow = 2,
  rel_heights = c(1.25,1)
)

# Height follows the submitted Fig 1 aspect (h/w ~= 1.466, taller than the 9.5 in
# page cap would give); width is the Science 2-column 7.25 in.
fig1_h <- round(page_w * 1.466, 2)   # ~10.63 in

save_figure_3way(fig_1, name = "figure1",
                 width = page_w, height = fig1_h,
                 dir = file.path(here::here(), "figures"))
