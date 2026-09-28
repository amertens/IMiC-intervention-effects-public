# =============================================================================
# figS2-milq-deficiency.R
#
# Builds Fig S2: relative risk (intervention vs control, 95% CI, log scale) of a milk
# nutrient being deficient, i.e. below the 10th percentile of the age-specific MILQ
# reference, one facet per study x visit. Points are coloured by significance
# (grey = not significant, blue = P < 0.05 before FDR, orange = FDR-significant).
#
# Inputs:  results/milq_deficiency_reduction_analysis_results.RDS
#            (src/2 analysis/milq-deficiency-reduction-analysis.R)
# Outputs: figures/figureS2_milq_deficiency.png
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

# canonical biomarker labels (0_figure-functions.R)
plot_df$label_f <- canonical_label(plot_df$biomarker, plot_df$label_f)

plot_df <- plot_df %>% mutate(
  studytime=as.character(studytime),
  studytime = case_when(
    studytime=="Misame-1" ~ "MISAME-III (14-21 days)",
    studytime=="Misame-2" ~ "MISAME-III (1-2 mo.)",
    studytime=="Misame-3" ~ "MISAME-III (3-4 mo.)",
    studytime=="Vital-40" ~ "Mumta-LW (1.5 mo.)",
    studytime=="Vital-56" ~ "Mumta-LW (2 mo.)",
    studytime=="Elicit-1" ~ "ELICIT (1 mo.)",
    studytime=="Elicit-5" ~ "ELICIT (5 mo.)"
  ),
  study = case_when(
    grepl("MISAME", studytime) ~ "MISAME-III",
    grepl("Mumta", studytime) ~ "Mumta-LW",
    grepl("ELICIT", studytime) ~ "ELICIT"
  ),
  studytime=factor(studytime, levels=c("MISAME-III (14-21 days)", "MISAME-III (1-2 mo.)", "MISAME-III (3-4 mo.)",
                                       "Mumta-LW (1.5 mo.)", "Mumta-LW (2 mo.)", "ELICIT (1 mo.)", "ELICIT (5 mo.)"))
)

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
  "ELICIT (1 mo.)",
  "ELICIT (5 mo.)",
  "Mumta-LW (1.5 mo.)",
  "Mumta-LW (2 mo.)",
  "MISAME-III (14-21 days)",
  "MISAME-III (1-2 mo.)",
  "MISAME-III (3-4 mo.)"
)

plot_df <- plot_df %>%
  mutate(
    studytime = factor(studytime, levels = studytime_levels)
  )

SIG_KEY <- c("Not Significant" = "Not significant", "Significant before FDR" = "Nominally significant",
             "Significant after FDR" = "FDR-significant")
p <- ggplot(plot_df, aes(x=label_f , y=RR, color=factor(sigcat), alpha=factor(sigcat), shape=factor(sigcat))) +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype="dashed") +
  facet_wrap(~studytime, scale="free_x", nrow=2,  drop = FALSE) +
  # every other doubling: at eight breaks "1/16" and "1/8" overprinted as "1/161/8"
  scale_y_log10(  breaks = c(1/16, 1/4, 1, 4),
                  labels = c("1/16", "1/4", "1", "4")) +
  coord_flip(    ylim = c(0.025, 8),
                 clip = "on") +
  ggtitle("") +
  theme_imic(base_size = 8) +   # shared theme (Helvetica, font-size floors)
  # Named mappings so each significance level keeps its colour/alpha/shape regardless of
  # which levels are present in the data (an unnamed scale shifts "Significant after FDR"
  # onto blue when "Significant before FDR" is absent). Sig-after-FDR is always orange.
  # Shape matches Fig 2 / Fig S3's 3-tier scheme: open circle for not-significant and
  # significant-before-FDR, filled circle only for FDR-significant. Key worded as in Fig 2.
  scale_alpha_manual(values=c("Not Significant"=0.3, "Significant before FDR"=0.5, "Significant after FDR"=0.8),
                     labels = SIG_KEY, name = "Statistical Significance") +
  scale_color_manual(values=c("Not Significant"="grey50", "Significant before FDR"=tableau10[1], "Significant after FDR"=tableau10[2]),
                     labels = SIG_KEY, name = "Statistical Significance") +
  scale_shape_manual(values=c("Not Significant"=1, "Significant before FDR"=1, "Significant after FDR"=19),
                     labels = SIG_KEY, name = "Statistical Significance") +
  # theme_imic()'s white strips, with plain (unbolded) titles as in the main figures
  theme(legend.position = "inside",
        legend.position.inside = c(0.87,0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "plain"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3)) +
  ylab("Relative Risk") + xlab("Milk component")

# ragg::agg_png renders UTF-8 glyphs (α/γ tocopherol, B-vitamin subscripts); the default
# Windows png device cannot (mbcsToSbcs conversion failure).
ggsave(here("figures/figureS2_milq_deficiency.png"), p, width = 8, height = 6, dpi = 300, device = ragg::agg_png)


