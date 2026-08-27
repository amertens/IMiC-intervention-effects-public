# =============================================================================
# figure-scripts/manuscript_figures/figS3-milq-deficiency.R
#
# Reads:  results/milq_deficiency_reduction_analysis_results.RDS
# Writes: figure-data/fig-milq-deficiency-reduction-data.RDS
#         figure-data/fig-milq-deficiency-reduction-data.xlsx
#         figures/fig-milq-deficiency-reduction-forest-plot.png
#
# Paths above were recovered from this script's syntax tree and are
# repo-relative; they resolve from the repo root via here::here().
#
# Header generated from the code itself; it makes no claim about method.
# See README.md for run order and results/ARTIFACT_MANIFEST.csv for the
# exhibit each script feeds.
# =============================================================================



rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))


#load results
df = readRDS(paste0(here::here(),"/results/milq_deficiency_reduction_analysis_results.RDS"))

df= clean_biomarker_labels(df)

plot_df <- df %>%
  mutate(sigcat=case_when(
    sigFDR==1 ~ "Significant after FDR",
    sigFDR==0 & sig ==1 ~ "Significant before FDR",
    sigFDR==0 & sig ==0 ~ "Not Significant"
  ),
  sigcat=factor(sigcat, levels=rev(c("Significant after FDR","Significant before FDR","Not Significant"))))

temp =plot_df %>% filter(biomarker=="se") %>% mutate(est_f=paste0(studytime,": ",round(RR,3)," (",round(ci.lb,3),", ",round(ci.ub,3),")"))
temp$est_f

# canonical biomarker labels (single source of truth: 0_figure-functions.R)
plot_df$label_f <- canonical_label(plot_df$biomarker, plot_df$label_f)

plot_df <- plot_df %>% mutate(
  studytime=as.character(studytime),
  studytime = case_when(
    studytime=="Misame-1" ~ "MISAME-III\n(14-21 days)",
    studytime=="Misame-2" ~ "MISAME-III\n(1-2 mo.)",
    studytime=="Misame-3" ~ "MISAME-III\n(3-4 mo.)",
    studytime=="Vital-40" ~ "Mumta-LW\n(1.5 mo.)",
    studytime=="Vital-56" ~ "Mumta-LW\n(2 mo.)",
    studytime=="Elicit-1" ~ "ELICIT\n(1 mo.)",
    studytime=="Elicit-5" ~ "ELICIT\n(5 mo.)"
  ),
  study = case_when(
    grepl("MISAME", studytime) ~ "MISAME-III",
    grepl("Mumta", studytime) ~ "Mumta-LW",
    grepl("ELICIT", studytime) ~ "ELICIT"
  ),
  studytime=factor(studytime, levels=c("MISAME-III\n(14-21 days)", "MISAME-III\n(1-2 mo.)", "MISAME-III\n(3-4 mo.)",
                                       "Mumta-LW\n(1.5 mo.)", "Mumta-LW\n(2 mo.)", "ELICIT\n(1 mo.)", "ELICIT\n(5 mo.)"))
)

dput(unique(plot_df$label_f))
plot_df <- plot_df %>% 
  mutate(label_f=factor(label_f, levels=rev(c("Total Energy",  "Total Fat", "Total Protein",
                                          "Vitamin B₂", "Vitamin B₃",  "Vitamin B₅",
                                          "Vitamin B₆", "Vitamin B₇","Vitamin B₁₂",
                                          "Vitamin A","α-Tocopherol", "γ-Tocopherol",
                                          "Calcium", "Choline", "Copper","Iron",
                                          "Magnesium", "Phosphorus", "Potassium",
                                          "Selenium", "Sodium",
                                          "Zinc"  ))))

studytime_levels <- c(
  "ELICIT\n(1 mo.)",
  "ELICIT\n(5 mo.)",
  "Mumta-LW\n(1.5 mo.)",
  "Mumta-LW\n(2 mo.)",
  "MISAME-III\n(14-21 days)",
  "MISAME-III\n(1-2 mo.)",
  "MISAME-III\n(3-4 mo.)"
)

plot_df <- plot_df %>%
  mutate(
    studytime = factor(studytime, levels = studytime_levels)
  )

p <- ggplot(plot_df, aes(x=label_f , y=RR, color=factor(sigcat), alpha=factor(sigcat), shape=factor(sigcat))) +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype="dashed") +
  facet_wrap(~studytime, scale="free_x", nrow=2,  drop = FALSE) +
  #facet_grid(~studytime) +
  scale_y_log10(  breaks = c(0.125/2,0.125, 0.25, 0.5, 1, 2, 4, 8),
                  labels = c("1/16","1/8","1/4", "1/2", "1", "2", "4", "8")) +
  coord_flip(    ylim = c(0.025, 8),
                 clip = "on") +
  ggtitle("") +
  theme_imic(base_size = 8) +   # Science-submission theme (Helvetica, font floors); harmonized 2026-08-26
  # NAMED mappings so each significance level keeps its colour/alpha/shape regardless of
  # which levels are present in the data (an unnamed scale shifts "Significant after FDR"
  # onto blue when "Significant before FDR" is absent). Sig-after-FDR is always orange.
  # Shape matches Fig 2 / Fig S3's 3-tier scheme: open circle for not-significant and
  # significant-before-FDR, filled circle only for FDR-significant (harmonized 2026-08-26;
  # previously this panel used solid circles for every tier).
  scale_alpha_manual(values=c("Not Significant"=0.3, "Significant before FDR"=0.5, "Significant after FDR"=0.8)) +
  scale_color_manual(values=c("Not Significant"="grey50", "Significant before FDR"=tableau10[1], "Significant after FDR"=tableau10[2])) +
  scale_shape_manual(values=c("Not Significant"=1, "Significant before FDR"=1, "Significant after FDR"=19)) +
  # strip.background/strip.text overrides removed (harmonized 2026-08-26): theme_imic's
  # own bold grey90 strip now applies, matching the other harmonized figures.
  theme(legend.position = "inside",
        legend.position.inside = c(0.87,0.2),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3)) +
  ylab("Relative Risk") + xlab("Milk component")
p

# ragg::agg_png renders UTF-8 glyphs (α/γ tocopherol, B-vitamin subscripts); the default
# Windows png device cannot (mbcsToSbcs conversion failure).
ggsave(here("figures/fig-milq-deficiency-reduction-forest-plot.png"), p, width = 8, height = 6, dpi = 300, device = ragg::agg_png)
saveRDS(plot_df, file= paste0(here::here(),"/figure-data/fig-milq-deficiency-reduction-data.RDS"))
writexl::write_xlsx(plot_df, paste0(here::here(),"/figure-data/fig-milq-deficiency-reduction-data.xlsx"))


