# =============================================================================
# fig4-milq-boxplots.R
#
# Builds Fig 4: for 15 primary milk nutrients, violin plots of the participant-level
# concentrations (native units, log y-axis) by study x visit and arm (control vs
# pooled intervention), with the age-specific MILQ reference P50 (grey) and P10 (red)
# drawn over each visit. Asterisks mark combined-arm ATEs that are FDR-significant
# (* Q < 0.05, ** < 0.01, *** < 0.001).
#
# Inputs:  data/merged_analysis_datasets.RDS (participant-level)
#          data/milq_age_specific_cutoffs_clean.RDS (MILQ reference percentiles)
#          metadata/milk_component.Rdata
#          results/adjusted_combined_arms_intervention_effects_results_clean.RDS
# Outputs: figures/figure4.jpeg
# [needs restricted data] merged_analysis_datasets.RDS is participant-level data.
# =============================================================================

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

#--------------------------------------------------------------------------------
# load data
#--------------------------------------------------------------------------------

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
dfull <- readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS")) %>%
  rename_with(~ "vitamin.a", any_of("vitamin.A"))   # normalize FSV vitamin A column casing (robust to either)
int_res <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))

# MILQ age-specific reference percentiles (P10, P50), converted to the IMiC units:
# x1000 for B3, Cu, Fe, B5 and Zn; /10 for protein, fat and carbohydrate; B12 from
# ug/L to pmol/L as x 1e6 / 1355 (B12 molar mass 1355 g/mol).
milq <- readRDS(here("data/milq_age_specific_cutoffs_clean.RDS"))
milq <- milq %>%
  select(nutrient, age_band, P10, P50) %>%
  mutate(P10 = case_when(
    nutrient=="b3" ~ P10*1000,
    nutrient=="b12" ~ P10 * 1e6 / 1355,   # ug/L -> pmol/L
    nutrient=="cu" ~ P10*1000,
    nutrient=="fe" ~ P10*1000,
    nutrient=="pa" ~ P10*1000,
    nutrient=="zn" ~ P10*1000,
    nutrient=="protein" ~ P10/10,
    nutrient=="fat" ~ P10/10,
    nutrient=="carbohydrate" ~ P10/10,
    nutrient==nutrient ~ P10),
    P50 = case_when(
      nutrient=="b3" ~ P50*1000,
      nutrient=="b12" ~ P50 * 1e6 / 1355,   # ug/L -> pmol/L
      nutrient=="cu" ~ P50*1000,
      nutrient=="fe" ~ P50*1000,
      nutrient=="pa" ~ P50*1000,
      nutrient=="zn" ~ P50*1000,
      nutrient=="protein" ~ P50/10,
      nutrient=="fat" ~ P50/10,
      nutrient=="carbohydrate" ~ P50/10,
      nutrient==nutrient ~ P50)) %>%
  dplyr::rename(biomarker=nutrient) %>%
  # lower-case to match the pivot_longer'd `d$biomarker` (and int_res$biomarker,
  # already lower-cased below); "vitamin.A" vs "vitamin.a" silently failed this
  # join and dropped every Vitamin A row via the downstream filter(!is.na(P10)).
  mutate(biomarker = str_to_lower(biomarker))




int_res <- int_res %>% filter(measure=="ATE") %>%
  # Use the contrast p-value (pval) and its FDR-adjusted q (pval_adj). The
  # chi_pval / chi_pval_adj columns are the omnibus chi-square test of any intervention
  # effect and are not appropriate for these significance markers.
  select(study, visit,  chi_pval, biomarker, chi_pval_adj, pval, pval_adj, sig, sigFDR ) %>%
  mutate(  visit=case_when(
    study=="Elicit" & visit=="1 mo."   ~ "1 month",
    study=="Elicit" & visit=="5 mo."   ~ "5 months",
    study=="Vital" & visit=="1.5 mo."   ~ "1.5 months",
    study=="Vital" & visit=="2 mo."   ~ "2 months",
    study=="Misame" & visit=="14-21 days"  ~ "14-21 days",
    study=="Misame" & visit=="1-2 mo." ~ "1-2 months",
    study=="Misame" & visit=="3-4 mo." ~ "3-4 months"))

dfull <- dfull %>% mutate(
  arm_f=arm,
  visit=case_when(
    studyid=="ELICIT" & visit=="1"   ~ "1 month",
    studyid=="ELICIT" & visit=="5"   ~ "5 months",
    studyid=="VITAL-Lactation" & visit=="40"   ~ "1.5 months",
    studyid=="VITAL-Lactation" & visit=="56"   ~ "2 months",
    studyid=="MISAME-3" & visit=="1"   ~ "14-21 days",
    studyid=="MISAME-3" & visit=="2"   ~ "1-2 months",
    studyid=="MISAME-3" & visit=="3"   ~ "3-4 months"),
  arm = case_when(
    arm=="Az." ~ "Control",
    arm=="BEP+ExBf+AZT" ~ "Intervention",
    arm=="BEP+ExBf" ~ "Intervention",
    arm=="Nico" ~ "Intervention",
    arm=="Nico+Az." ~ "Intervention",
    arm=="BEP/BEP" ~ "Intervention",
    arm=="IFA/BEP" ~ "Intervention",
    arm=="BEP/IFA" ~ "Control",
    arm==arm ~ arm),
  arm_visit=paste0(arm, " ", visit)) %>%
  mutate(arm=factor(arm, levels=c("Control","Intervention")),
         arm_visit=factor(arm_visit,
                          levels=rev(c( "Control 14-21 days","Intervention 14-21 days", "Control 1-2 months", "Intervention 1-2 months",  "Control 3-4 months",
                                        "Intervention 3-4 months",
                                        "Control 1 month", "Intervention 1 month", "Control 5 months", "Intervention 5 months",
                                        "Control 1.5 months","Intervention 1.5 months", "Control 2 months",  "Intervention 2 months"))))


# MILQ age bands matched to each study's visits
milq  <- bind_rows(milq  %>% mutate(studyid="MISAME-3"),
                   milq  %>% mutate(studyid="ELICIT"),
                   milq  %>% mutate(studyid="VITAL-Lactation")) %>%
  mutate(
    visit= case_when(
      studyid=="MISAME-3" & age_band=="4-17 d" ~ "14-21 days",
      studyid=="MISAME-3" & age_band=="1-2 m" ~ "1-2 months",
      studyid=="MISAME-3" & age_band=="3-4 m" ~ "3-4 months",
      studyid=="ELICIT" & age_band=="1-2 m" ~ "1 month",
      studyid=="ELICIT" & age_band=="4-5 m" ~ "5 months",
      studyid=="VITAL-Lactation" & age_band=="1-2 m" ~ "1.5 months",
      studyid=="VITAL-Lactation" & age_band=="2-3 m" ~ "2 months")
  ) %>% filter(!is.na(visit))

d <- dfull %>% select(study, studyid, visit, arm, arm_visit,
                      b3, nmn, nufa, nr, nam, nad, pa, ribo, t,
                      b1, fmn, tmp, b2, pl, pm, b6, se, g.tocopherol,
                      cu, trp, b12, fad, bio, a.tocopherol, vitamin.a,
                      kcal.l, fat, protein) %>%
  #transform to long
  pivot_longer(cols = c(b3, nmn, nufa, nr, nam, nad, pa, ribo, t,
                        b1, fmn, tmp, b2, pl, pm, b6, se, g.tocopherol,
                        cu, trp, b12, fad, bio, a.tocopherol, vitamin.a,
                        kcal.l, fat, protein), names_to = "biomarker", values_to = "value")

# keep only nutrients with a MILQ reference
d <- left_join(d, milq, by=c("studyid", "biomarker", "visit")) %>% filter(!is.na(P10))

d$visit <- gsub(" months","", d$visit)
d$visit <- gsub(" month","", d$visit)
d$visit <- gsub(" days","", d$visit)
int_res$visit <- gsub(" months","", int_res$visit)
int_res$visit <- gsub(" month","", int_res$visit)
int_res$visit <- gsub(" days","", int_res$visit)

int_res$biomarker <- str_to_lower(int_res$biomarker)

# lower-case the join key so the merged-data column "vitamin.A" matches the cleaned-results
# "vitamin.a" (case mismatch otherwise drops vitamin A's significance to NA).
df <- left_join(d %>% mutate(biomarker = str_to_lower(biomarker)), int_res, by=c("study","visit","biomarker"))


#clean up labels for the plot
# curated set of biomarkers shown in Fig 4 (panel order matches yaxis_lab_vec below);
# display labels come from the canonical map in 0_figure-functions.R so they match
# every other figure (α/γ-tocopherol, B-vitamin subscripts, etc.).
displayed_codes <- c("kcal.l","fat","protein","b1","b2","b3","pa","b6","bio","b12",
                     "vitamin.a","a.tocopherol","g.tocopherol","se","cu")
d <- df %>% mutate(study_f =
                     case_when(
                       study == "Elicit" ~ "El",
                       study == "Vital" ~ "Mu",
                       study == "Misame" ~ "Mi"
                     ),
                   visit_f =
                     case_when(
                       visit == "14-21" ~ "0.5",
                       visit == "1-2" ~ "1.5",
                       visit == "3-4" ~ "3.5",
                       visit == "1.5" ~ "1.5",
                       visit == "2" ~ "2",
                       visit == "1" ~ "1",
                       visit == "5" ~ "5"
                     ),
                   study_f = factor(study_f, levels = c("Mu","Mi",  "El")),
                   visit_f = factor(visit_f, levels = c("0.5", "1.5", "2", "3.5", "1", "5")),
                   study_visit=paste0(study_f,"-", visit_f),
                     biomarker_f = ifelse(biomarker %in% displayed_codes,
                                          canonical_label(biomarker), NA_character_)
                   )


# 0. Helper: asterisks from p-value  --------------------------------------------
getAsterisks <- function(p) {
  if (is.na(p))         return("")
  if (p < 0.001)        return("***")
  if (p < 0.01)         return("**")
  if (p < 0.05)         return("*")
  ""
}

# 1.  Add the asterisk column once  --------------------------------------------
# asterisks reflect the FDR-adjusted q-value of the contrast (pval_adj), not the omnibus chi_pval.
d <- d %>% mutate(asterisks = map_chr(pval_adj, getAsterisks))

# 2.  Colours for the two arms --------------------------------------------------
# the blue/orange of Fig 2 (tableau10 1-2)
arm_cols <- c("Control" = "#1F77B4",
              "Intervention" = "#FF7F0E")


# One panel per nutrient: violins by study x visit and arm, MILQ P50/P10 segments,
# significance asterisks above each visit.
plot_biomarker <- function(biomarker_name,
                           d,
                           arm_cols,
                           yaxis_lab="",
                           theme_imic) {

  ## subset data for biomarker
  dat <- d %>% filter(biomarker_f == biomarker_name)

  ## create label_df with significance stars
  lab_df <- dat %>%
    group_by(study, visit, study_visit, biomarker) %>%
    summarise(
      x_lab = first(study_visit),
      y_lab = fivenum(value, na.rm = TRUE)[5] * 1,
      label = first(asterisks),
      .groups = "drop"
    ) %>%
    filter(label != "") %>% distinct()

  ## study-specific axis tick colours, the study colours of study_colors.R. The
  ## study code is part of each tick label, so colour here is a redundant cue.
  study_colors <- c(
    "Mi" = "#E69F00",  # MISAME-III (Okabe-Ito orange)
    "Mu" = "#CC79A7",  # Mumta-LW / VITAL (Okabe-Ito reddish purple)
    "El" = "#0072B2"   # ELICIT (Okabe-Ito blue)
  )

  ## fixed ordering of study_visit
  dat$study_visit <- factor(
    dat$study_visit,
    levels = c(
      "El-1", "El-5",
      "Mi-0.5", "Mi-1.5", "Mi-3.5",
      "Mu-1.5", "Mu-2"
    )
  )

  xlabs <- levels(dat$study_visit)
  xcols <- study_colors[substr(xlabs, 1, 2)]
  # "El-1" -> El with subscript 1 (no hyphen)
  sub_lab <- function(x) parse(text = sprintf('"%s"["%s"]', substr(x, 1, 2), substring(x, 4)),
                               keep.source = FALSE)

  ## MILQ P50 / P10 reference line segments
  ref_df <- dat %>%
    group_by(study_visit) %>%
    summarise(P10 = unique(P10),
              P50 = unique(P50),.groups = "drop") %>%
    mutate(
      x_start = as.numeric(study_visit) - 0.4,
      x_end   = as.numeric(study_visit) + 0.4
    )

  ## main plot
  p <- ggplot(dat, aes(x = study_visit, y = value, fill = arm)) +
    # no outline (colour = NA), black median line only
    geom_violin(colour = NA,
                quantiles = 0.5, quantile.linetype = "solid",
                quantile.colour = "black", quantile.linewidth = 0.4,
                trim = TRUE,
                position = position_dodge(width = 0.8)) +

    geom_segment(
      data = ref_df,
      aes(x = x_start, xend = x_end, y = P50, yend = P50),
      inherit.aes = FALSE,
      linewidth = 0.75,
      colour = "#7F7F7F"
    ) +
    geom_segment(
      data = ref_df,
      aes(x = x_start, xend = x_end, y = P10, yend = P10),
      inherit.aes = FALSE,
      linewidth = 0.75,   # same weight as the grey P50 line
      colour = "#D62728"
    ) +

    scale_y_log10() +
    scale_x_discrete(labels = sub_lab) +
    scale_fill_manual(values = arm_cols) +
    scale_colour_manual(values = arm_cols) +
    theme_imic() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(colour = xcols),
      plot.title = element_text(face = "plain", hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3)
    ) +
    labs(
      title = biomarker_name,
      y = yaxis_lab,
      fill = "Arm"
    ) +
    geom_text(
      data = lab_df,
      aes(x = x_lab, y = y_lab, label = label),
      inherit.aes = FALSE,
      size = 4, vjust = 1
    )
  return(p)
}


# canonical labels, same order as displayed_codes / yaxis_lab_vec
biomarkers_to_plot <- canonical_label(displayed_codes)

yaxis_lab_vec <- c(
  "kcal/L",     # Energy density
  "g/100 mL",   # Fat (spec unit "percent")
  "g/100 mL",   # Protein (spec unit "percent")

  "\u00b5g/L",       # Vitamin B1
  "\u00b5g/L",       # Vitamin B2
  "\u00b5g/L",  # Vitamin B3 (as nicotinamide; spec unit ug/L)
  "\u00b5g/L",  # Vitamin B5 (pantothenic acid; spec unit ug/L)
  "\u00b5g/L",       # Vitamin B6
  "\u00b5g/L",       # Vitamin B7 (biotin)
  "pmol/L",     # Vitamin B12 (IMiC scale)

  "mg/L",       # Vitamin A (retinol)
  "mg/L",       # Alpha-tocopherol
  "mg/L",       # Gamma-tocopherol

  "\u00b5g/L",       # Selenium
  "\u00b5g/L"   # Copper (median ~430 ug/L)
)

plot_list <- Map(
  function(biomarker, ylab) {
    plot_biomarker(
      biomarker_name = biomarker,
      d = d,
      yaxis_lab = ylab,
      arm_cols = arm_cols,
      theme_imic = theme_imic
    )
  },
  biomarkers_to_plot,
  yaxis_lab_vec
)

library(gridExtra)

fig4 <-grid.arrange(
  grobs = plot_list,
  ncol = 3
)
# ragg::agg_jpeg renders the UTF-8 panel titles (α/γ-tocopherol, B-vitamin subscripts).
ggsave(fig4,
       filename = here("figures/figure4.jpeg"),
       width = 10, height = 8, units = "in", dpi = 300, device = ragg::agg_jpeg)
