# ---------------------------------------------------------------------------
# Auto-converted from the corresponding .Rmd (knitr::purl). Run via Rscript.
# Sets the working directory to figure-scripts/manuscript_figures/ so the
# script's '../../' relative paths resolve to the repo root.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(here))
setwd(file.path(here::here(), "figure-scripts/manuscript_figures"))


library(tidyverse)
library(cowplot)
library(magick)   
library(tinytex)
library(tidygraph)
library(ggraph)

# Science aesthetic: load shared palette (tableau10), theme_imic() floors, and
# save_figure_3way() (PDF + EPS + PNG export for the typesetter).
source(file.path(here::here(), "figure-scripts/0_figure-functions.R"))

knitr::opts_chunk$set(
  echo       = FALSE,
  warning    = FALSE,
  message    = FALSE,
  fig.width  = 8.27,         # should be exactly A4
  fig.height = 11.69,
  out.width  = "\\paperwidth",
  out.height = "\\paperheight",
  fig.align  = "center",
  fig.pos    = "!h"          
)

# ── automatic page break after any chunk with 'pagebreak=TRUE' ───────────
knitr::knit_hooks$set(
  pagebreak = function(before, options, envir) {
    if (!before) return("\\newpage")  # insert after the chunk runs
  }
)

# ── Science 2-column full-page dimensions (§3.11) ────────────────────────
page_w <- science_dims$full_page$width        # 7.25 in (184 mm, Science 2-col)
page_h <- science_dims$full_page$height_max   # 9.5  in (24 cm Science max page height)

#blank plots
blank_plot <- ggplot() + theme_void()

# ── helper: wrap a PNG as a ggdraw() canvas ─────────────────────────────
png_to_ggdraw <- function(path){
  ggdraw() + draw_image(image_read(path))
}


tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")





#Combined SL and PCA plots?

# Full SL variable-importance plot data, which INCLUDES the ELICIT infant-
# azithromycin NEGATIVE-CONTROL arm (built by src/3 visualizations/5-SL_VIM_plots.R).
SLvim_all <- readRDS("../../figure-data/SL_vim_plot_data.RDS")

# Shared CV-AUC x-range that covers EVERY confidence interval across the three
# studies (some ci.lb reach ~0.26), so nothing is clipped. Data-driven with
# minimal padding to keep white space down.
slvim_xlim <- c(floor(min(SLvim_all$ci.lb, na.rm = TRUE) / 0.05) * 0.05,
                min(1.0, ceiling(max(SLvim_all$ci.ub, na.rm = TRUE) / 0.05) * 0.05))

# ---- Rebuild Panel A from underlying data ---------------------------------
# The ggplot objects stored in SLvim_lab were built with an older ggplot2
# version and fail grid.draw under current ggplot2; rebuild from each plot's
# stored $data to guarantee a clean render.
build_slvim_panel <- function(d, title = "", show_x = FALSE,
                              legend_pos = "none", xlab_text = "",
                              xlim_range = slvim_xlim) {
  # Order predictor groups so "All" stays at the top of the y-axis
  group_levels <- c("All", "Macronutrients", "Micronutrients",
                    "B-vitamins", "HMOs", "Targeted proteins",
                    "Targeted metabolomics")
  d$group <- factor(d$group, levels = rev(intersect(group_levels, unique(d$group))))
  ggplot(d, aes(x = cvAUC, y = group, colour = visit_f)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = ci.lb, xmax = ci.ub),
                   position = position_dodge(width = 0.6)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
    facet_grid(. ~ arm_f) +
    scale_colour_manual(
      values = c("<1 month" = "#7F7F7F",     # grey / orange / teal, matching the submission
                 "1-2 months" = "#FF7F0E",
                 "2-5 months" = "#17BECF"),
      name = "Collection time") +
    scale_x_continuous(breaks = c(0.4, 0.6, 0.8, 1.0),
                       expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(xlim = xlim_range) +   # zoom, not filter -> keeps every CI
    labs(x = if (show_x) "CV-AUC" else xlab_text, y = title) +
    # base_family added for Helvetica consistency (2026-08-26). Strip background kept
    # WHITE/blank (reverted from a gray90 fill tried the same day) per author preference
    # -- Fig 1 should match Fig 2's white-strip forest style, not theme_imic()'s gray90.
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 8),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          legend.position  = legend_pos,
          legend.text      = element_text(size = 7),
          legend.title     = element_text(size = 7))
}
# Authoritative PCA-of-modality ATE results (all 8 milk modalities, INCLUDING
# Untargeted metabolomics + Microbiome). The old figures/figure-data/ copy was a
# stale 6-modality subset (dropped untargeted + microbiome), so Panel B was missing
# two rows; src/2 analysis/3_adjusted_analysis_pca.R writes the full 8 here.
pca_df <- readRDS(file="../../results/pca_intervention_effects_results.RDS") %>% filter(measure=="ATE")
pca_df$studytime <- gsub("-40", "-1.5", pca_df$studytime)
pca_df$studytime <- gsub("-56", "-2", pca_df$studytime)
pca_df$studytime <- gsub("Vital", "Mumta-LW", pca_df$studytime)
pca_df$studytime <- gsub("Misame-1", "Misame (14-21 days)", pca_df$studytime)
pca_df$studytime <- gsub("Misame-2", "Misame (1-2 mo)", pca_df$studytime)
pca_df$studytime <- gsub("Misame-3", "Misame (3-4 mo)", pca_df$studytime)
pca_df$studytime <- gsub("-1.5", " (1.5 mo.)", pca_df$studytime)
pca_df$studytime <- gsub("Elicit-1", "ELICIT (1 mo.)", pca_df$studytime)
pca_df$studytime <- gsub("Mumta-LW-2", "Mumta-LW (2 mo.)", pca_df$studytime)
pca_df$studytime <- gsub("Elicit-5", "ELICIT (5 mo.)", pca_df$studytime)

pca_df$studytime



p_elicit <- build_slvim_panel(dplyr::filter(SLvim_all, studyid == "ELICIT"),    title = "ELICIT")
p_vital  <- build_slvim_panel(dplyr::filter(SLvim_all, studyid == "Mumta-LW"),  title = "Mumta-LW")
p_misame <- build_slvim_panel(dplyr::filter(SLvim_all, studyid == "Misame-III"), title = "MISAME-III",
                              show_x = TRUE, legend_pos = "bottom")

SLvim_lab_plot <- plot_grid(p_elicit, p_vital, p_misame,
                            ncol = 1, rel_heights = c(1, 1, 1.25),
                            align = "v", axis = "lr")
SLvim_lab_plot




# p_pca <- ggplot(pca_df, aes(x=label_f, y=est, color=label_f)) + geom_point() +
#   geom_linerange(aes(ymin=cil, ymax=ciu)) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_flip() +
#   facet_wrap(studytime~contrast, scale="free") +
#   #ggtitle("Intervention Effects on the First Principal Component") +
#   scale_color_manual(values=rev(tableau10)) +
#   #theme_imic() +
#   theme(strip.background = element_blank(),
#         axis.text = element_text(size = 7),
#         strip.text = element_text(size = 8),
#         legend.position="none") + 
#   ylab("Average Treatment Effect") + xlab("Milk modality")
# 
# p_pca

pca_df <- pca_df %>% mutate(
  studytime=gsub("Misame","MISAME-III",studytime),
  studyid=case_when(
    grepl("ELICIT",studytime) ~"ELICIT",
    grepl("Mumta-LW",studytime) ~"Mumta-LW",
    grepl("MISAME",studytime) ~"MISAME-III"
  )#,
  # studytime=case_when(
  #   grepl("ELICIT",studytime) ~paste0(studytime,"\nNicotinamide"),
  #   grepl("Mumta-LW",studytime) ~paste0(studytime,"\nBEP"),
  #   grepl("MISAME",studytime) ~paste0(studytime,"\nBEP"),
  # )
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
  geom_hline(yintercept = 0, linetype="dashed") +
  coord_flip() +
  facet_grid(~studytime) +
  scale_color_manual(values=rev(tableau10)) +
  # base_family added for Helvetica consistency; strip stays WHITE/blank (see
  # build_slvim_panel() above -- reverted from gray90 same day, per author preference).
  theme_bw(base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 7),
        strip.text = element_text(size = 8),
        legend.position="none") +
  ylab("") + xlab("")


p_pca2 <- ggplot(pca_df %>% filter(studyid=="Mumta-LW"), aes(x=label_f, y=est, color=label_f)) + geom_point() +
  geom_linerange(aes(ymin=cil, ymax=ciu)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  coord_flip() +
  facet_grid(~studytime) +
  scale_color_manual(values=rev(tableau10)) +
  # base_family added for Helvetica consistency; strip stays WHITE/blank (see
  # build_slvim_panel() above -- reverted from gray90 same day, per author preference).
  theme_bw(base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 7),
        strip.text = element_text(size = 8),
        legend.position="none") +
  ylab("") + xlab("Milk modality")

p_pca3 <- ggplot(pca_df %>% filter(studyid=="MISAME-III") %>%
                   mutate(studytime=factor(studytime, levels = c("MISAME-III (14-21 days)","MISAME-III (1-2 mo)","MISAME-III (3-4 mo)" ))), aes(x=label_f, y=est, color=label_f)) + geom_point() +
  geom_linerange(aes(ymin=cil, ymax=ciu)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  coord_flip() +
  facet_grid(~studytime) +
  scale_color_manual(values=rev(tableau10)) +
  # base_family added for Helvetica consistency; strip stays WHITE/blank (see
  # build_slvim_panel() above -- reverted from gray90 same day, per author preference).
  theme_bw(base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 7),
        strip.text = element_text(size = 8),
        legend.position="none") +
  ylab("Average Treatment Effect") + xlab("")


p_pca <- plot_grid(p_pca1 ,
                            p_pca2 , 
                            p_pca3,
                            ncol=1, rel_heights=c(1,1,1.1))
p_pca




fig_2 <- plot_grid(
  SLvim_lab_plot,
  p_pca,
  align="h",
  axis ="l",
      labels="AUTO",
  ncol = 1, nrow = 2,
  rel_heights = c(1.25,1)  
)
fig_2

# Taller height to match the submitted Fig 1 aspect (h/w ~= 1.466, vs 1.31 at the
# 9.5 in cap); width stays at the Science 2-column 7.25 in.
fig1_h <- round(page_w * 1.466, 2)   # ~10.63 in

#save as png (legacy) + 3-way (PDF/EPS/PNG) for Science submission
ggsave(
  filename = "../../figures/figure2_PCA.png",
  plot = fig_2,
  width = page_w, height = fig1_h,
  units = "in", dpi = 300
)
save_figure_3way(fig_2, name = "figure1",
                 width = page_w, height = fig1_h,
                 dir = file.path(here::here(), "figures"))


