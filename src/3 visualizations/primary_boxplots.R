# =============================================================================
# src/3 visualizations/primary_boxplots.R
#
# Reads:  data/merged_analysis_datasets.RDS
#         data/milq_age_specific_cutoffs_clean.RDS
#         metadata/milk_component.Rdata
#         results/adjusted_combined_arms_intervention_effects_results_clean.RDS
# Writes: figures/primary_boxplots.png
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


#NOTES! Fix merging (too many obs)
#and merge in MILQ cutoffs
#Make the 
#make the trim only affect the plotted points, not the boxplot

#--------------------------------------------------------------------------------
# plot themes
#--------------------------------------------------------------------------------

imic_palette  <- c(
  "grey90",      # non‑sig   (same as the qmd)
  tableau10[2]   # sig       (tableau10 is already defined in 0‑config.R)
)

theme_imic <- function(base_size = 8, base_family = "Helvetica") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.text.y        = element_text(size = base_size, hjust = 0),
      axis.text.x        = element_text(size = base_size),
      plot.title         = element_text(face = "bold"),
      strip.text         = element_text(face = "bold", size = base_size),
      legend.position    = "top",
      legend.title       = element_text(size = base_size),
      legend.text        = element_text(size = base_size * .9)
    )
}


#--------------------------------------------------------------------------------
# load data
#--------------------------------------------------------------------------------


load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
dfull <- readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))
int_res <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))
head(int_res)

# d[d$label_f=="nicotinamide riboside",]
# int_res%>% filter(biomarker=="Nr")
# 
# int_res <- int_res %>% filter( measure=="ATE", outcome_group  =="primary")
# sig_int_res <- int_res %>% filter(sigFDR==1)
# dput((unique(sig_int_res$biomarker)))

table(int_res$biomarker)

#B1 results to check:
res_check = int_res %>% filter(measure=="ATE",  biomarker=="B1") 


#Load Child milk component medians
#update to MILQ cutoffs
milq <- readRDS(here("data/milq_age_specific_cutoffs_clean.RDS"))
milq <- milq %>% filter(age_band=="1-2 m") %>%
  select(nutrient, P10) %>%
  mutate(P10 = case_when(
            nutrient=="b3" ~ P10*1000,
            nutrient=="cu" ~ P10*1000,
            nutrient=="fe" ~ P10*1000,
            nutrient=="pa" ~ P10*1000,
            nutrient=="zn" ~ P10*1000,
            nutrient=="protein" ~ P10/10,
            nutrient=="fat" ~ P10/10,
            nutrient=="carbohydrate" ~ P10/10)) %>%
  rename(biomarker=nutrient)

#convert units of variables





int_res <- int_res %>% filter(measure=="ATE") %>%
  select(study, visit,  chi_pval, biomarker, chi_pval_adj, sig, sigFDR ) %>%
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
table(dfull$arm_f, is.na(dfull$arm))
unique(dfull$arm_visit)


milq  <- bind_rows(milq  %>% mutate(studyid="MISAME-3"),
                   milq  %>% mutate(studyid="ELICIT"),
                   milq  %>% mutate(studyid="VITAL-Lactation"))




d <- dfull %>% select(study, studyid, visit, arm, arm_visit,
                      b3, nmn, nufa, nr, nam, nad, pa, ribo, t, 
                      b1, fmn, tmp, b2, pl, pm, b6, se, g.tocopherol, 
                      cu, trp, b12, fad, bio, a.tocopherol, vitamin.A, 
                      kcal.l, fat, protein) %>% 
  #transform to long
  pivot_longer(cols = c(b3, nmn, nufa, nr, nam, nad, pa, ribo, t, 
                        b1, fmn, tmp, b2, pl, pm, b6, se, g.tocopherol, 
                        cu, trp, b12, fad, bio, a.tocopherol, vitamin.A, 
                        kcal.l, fat, protein), names_to = "biomarker", values_to = "value") 

head(d)
head(milq)

dim(d)
dim(milq)
d <- left_join(d, milq, by=c("studyid", "biomarker"))
head(d)

table(d$biomarker, is.na(d$P10))

dim(d)
d$value

d <- d %>%
  filter(!is.na(P10)) %>%
  mutate(
    deficient =ifelse(value < P10, 1, 0)) 
table(d$deficient)
d<-clean_biomarker_labels(d)

head(d)


table(d$visit)
table(int_res$visit)

d$visit <- gsub(" months","", d$visit)
d$visit <- gsub(" month","", d$visit)
d$visit <- gsub(" days","", d$visit)
int_res$visit <- gsub(" months","", int_res$visit)
int_res$visit <- gsub(" month","", int_res$visit)
int_res$visit <- gsub(" days","", int_res$visit)

dim(d)
dim(int_res)
df <- left_join(d, int_res, by=c("study","visit","biomarker")) 
dim(df)




#clean up labels for the plot
d <- df %>% mutate(study_f =
                     case_when(
                       study == "Elicit" ~ "E",
                       study == "Vital" ~ "V",
                       study == "Misame" ~ "M"
                     ), 
                   study_f = factor(study_f, levels = c("M","V",  "E")),
                   visit = factor(visit, levels = c("14-21","1-2","3-4",  "1.5","2","1","5")))



# ──────────────────────────────────────────────────────────────
#  Multi‑vit figure: B1, B3, B12   (IMiC style, 10‑90 % trim)
# ──────────────────────────────────────────────────────────────
library(tidyverse)
library(ggpubr)     # geom_text() stars already present; ggpubr not strictly needed now
library(patchwork)  # wrap_plots()

# 0. Helper: asterisks from p‑value  --------------------------------------------
getAsterisks <- function(p) {
  if (is.na(p))         return("")
  if (p < 0.001)        return("***")
  if (p < 0.01)         return("**")
  if (p < 0.05)         return("*")
  ""
}

# 1.  Add the asterisk column once  --------------------------------------------
d <- d %>% mutate(asterisks = map_chr(chi_pval, getAsterisks))

# 2.  Trim to 10–90 % per biomarker/study/visit  -------------------------------
d_trim <- d %>% 
  group_by(biomarker, study, visit) %>%
  mutate(
    lo = quantile(value, 0.10, na.rm = TRUE),
    hi = quantile(value, 0.90, na.rm = TRUE),
    value_trim= ifelse(between(value, lo, hi), value, NA)) %>%
  ungroup()
summary(d_trim$value)
summary(d_trim$value_trim)

# 3.  Colours for the two arms --------------------------------------------------
arm_cols <- c("Control" = "#0072B2",
              "Intervention" = "#D55E00")



# 4.  Re‑usable plotting function ----------------------------------------------
outcome="B3"
plot_one_facet <- function(outcome) {
  
  # 4a. subset & label data
  dat <- d_trim %>% filter(biomarker == outcome)
  plot_label <- dat$label_f[1]
  
  lab_df <- dat %>% #filter() %>% # where to put the asterisks
    group_by(study, visit, biomarker) %>%     
    summarise(
      x_lab = interaction(study_f, visit, sep = " · ", drop = TRUE),
      y_lab = fivenum(value, na.rm = T)[4]  * 1.1,
      label = first(asterisks),
      .groups = "drop"
    ) %>% 
    filter(label != "") %>% distinct()
  
  # 4b. build ggplot ------------------------------------------------------------
  ggplot(dat,
         aes(x = interaction(study_f, visit, sep = " · "),
             y = value,
             fill = arm)) +
    geom_boxplot(outlier.shape = NA, alpha = .70,
                 position = position_dodge(width = .8)) +
    geom_jitter(aes(y=value_trim, colour = "black"),
                position = position_jitterdodge(jitter.width = .20,
                                                dodge.width  = .8),
                size = 1.2, alpha = .5, show.legend = FALSE) +
    geom_hline(aes(yintercept = P10), linetype = "dashed") +
    scale_y_log10() +
    scale_fill_manual(values = arm_cols) +
    scale_colour_manual(values = arm_cols) +
    theme_imic() +
    theme(legend.position = "none",
          axis.title.x = element_blank()) +
    labs(title = paste0(plot_label),
         y = "", fill = "Arm") +
    geom_text(data = lab_df,
              aes(x = x_lab, y = y_lab, label = label),
              inherit.aes = FALSE,
              size = 4, vjust = 0)
}

# 5.  graphical abstract plots 

plot_one_facet(outcome="B3")


# 5.  Generate the  plots --------------------------------------------------
sig_components <- c("B3", "Nmn", "Nufa", "Nr", "Nam", "Nad", "Pa", "Ribo", "T", 
                    "B1", "Fmn", "Tmp", "B2", "Pl", "Pm", "B6", "Se", "G.tocopherol", 
                    "Cu", "Trp", "B12", "Fad", "Bio", "A.tocopherol", "Vitamin.a", 
                    "Kcal.l", "Fat", "Protein")
plots <- map(sig_components, plot_one_facet)

# 6.  Arrange (2 columns) & print ----------------------------------------------
combined_plot <- wrap_plots(plots, ncol = 4)
print(combined_plot)


ggsave(combined_plot, file=here::here("figures/primary_boxplots.png"),
       width = 16, height = 12, units = "in")
