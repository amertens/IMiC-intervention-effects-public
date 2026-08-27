# =============================================================================
# figure-scripts/figures_s5_trajectory_lineplots.R
#
# Reads:  data/merged_analysis_datasets.RDS
#         metadata/milk_component.Rdata
#         results/combined_intervention_effects_results_combined_arms.RDS
# Writes: figure-data/figure_s4_trajectory_plots.RDS
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

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))
res <- readRDS(file=paste0(here::here(),"/results/combined_intervention_effects_results_combined_arms.RDS")) %>% filter(measure=="ATE")



d <- d %>%
  mutate(
    visit_f=case_when(
      visit==1 & study=="Misame" ~ "14-21 days",
      visit==2 & study=="Misame" ~ "1-2 mo.",
      visit==3 & study=="Misame" ~ "3-4 mo.",
      visit==1 & study=="Elicit" ~ "1 mo.",
      visit==5 & study=="Elicit" ~ "5 mo.",
      visit==40 & study=="Vital" ~ "1.5 mo.",
      visit==56 & study=="Vital" ~ "2 mo."
    )
  )
res <- res %>%
  mutate(
    visit_f=visit,
    visit=case_when(
      visit_f=="14-21 days" & study=="Misame" ~ 1,
      visit_f=="1-2 mo." & study=="Misame" ~ 2,
      visit_f=="3-4 mo." & study=="Misame" ~ 3,
      visit_f=="1 mo." & study=="Elicit" ~ 1,
      visit_f=="5 mo." & study=="Elicit" ~ 5,
      visit_f=="1.5 mo." & study=="Vital" ~ 40,
      visit_f=="2 mo." & study=="Vital" ~ 56
    )
  )


load(paste0(here::here(),"/metadata/milk_component.Rdata"))
d <- d %>% rename(vitamin.a=vitamin.A) %>%
  select(study, visit, arm, all_of(as.vector(unlist(all_milk_components)))) %>%
  #transform to long
  pivot_longer(cols=-c(study, visit, arm), names_to="biomarker", values_to="value") 
head(d)
head(res)

res$biomarker <- str_to_lower(res$biomarker)
d$visit
res$visit
table(d$visit)
table(res$visit)

dim(d)
dim(res)
d <- left_join(d, res, by=c("biomarker", "visit", "study")) 
dim(d)


#check missing categories
d$biomarker[is.na(d$category)]

d <- d %>% filter(!is.na(category))

unique(d$category)

#clean up categories
d$category[grepl("HMO",d$category)] <- "HMO" 
d$category[grepl("Amino",d$category)] <- "Amino Acid" 


head(d)
df <- d %>% group_by(study, visit,  biomarker) %>%
  mutate(
    value_raw=value,
    value =scale(value)
  ) %>% group_by(study, visit, arm, category) %>%
  summarise(
    value_raw =mean(value_raw, na.rm=T),
    value =mean(value, na.rm=T)
  ) %>%
  mutate(visit=factor(visit),
         arm=factor(arm))

#make lineplot with visit on x axis and value on y axis

# ggplot(df %>% filter(category=="B1"), aes(x=visit, y=value, group=arm, color=arm)) + 
#   geom_line() + geom_point() +
#   facet_wrap(~study, scales="free") +
#   
#   xlab("Visit") + ylab("Scaled value") +
#   ggtitle("Milk components by visit and category")


#to do:
# sort panels by number of sig effects
# make plots only averaging components with sig intervention effects

p_misame_scaled <- ggplot(df %>% filter(study=="Misame"), aes(x=visit, y=value, group=arm, color=arm)) + 
  geom_line() + geom_point() +
  facet_wrap(~category, scales="free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom") +
  xlab("Visit") + ylab("Scaled value") +
  ggtitle("A) Misame") 

p_elicit_scaled <- ggplot(df %>% filter(study=="Elicit"), aes(x=visit, y=value, group=arm, color=arm)) + 
  geom_line() + geom_point() +
  facet_wrap(~category, scales="free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom") +
  xlab("Visit") + ylab("Scaled value") +
  ggtitle("B) Elicit") 

p_vital_scaled <- ggplot(df %>% filter(study=="Vital"), aes(x=visit, y=value, group=arm, color=arm)) + 
  geom_line() + geom_point() +
  facet_wrap(~category, scales="free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom") +
  xlab("Visit") + ylab("Scaled value") +
  ggtitle("C) Vital") 



p_misame_unscaled <- ggplot(df %>% filter(study=="Misame"), aes(x=visit, y=value_raw, group=arm, color=arm)) + 
  geom_line() + geom_point() +
  facet_wrap(~category, scales="free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom") +
  xlab("Visit") + ylab("Value") +
  ggtitle("A) Misame") 

p_elicit_unscaled <- ggplot(df %>% filter(study=="Elicit"), aes(x=visit, y=value_raw, group=arm, color=arm)) + 
  geom_line() + geom_point() +
  facet_wrap(~category, scales="free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom") +
  xlab("Visit") + ylab("Value") +
  ggtitle("B) Elicit") 

p_vital_unscaled <- ggplot(df %>% filter(study=="Vital"), aes(x=visit, y=value_raw, group=arm, color=arm)) + 
  geom_line() + geom_point() +
  facet_wrap(~category, scales="free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom") +
  xlab("Visit") + ylab("Value") +
  ggtitle("C) Vital") 



saveRDS(
  list(p_misame_scaled=p_misame_scaled,
       p_elicit_scaled=p_elicit_scaled,
       p_vital_scaled=p_vital_scaled,
       p_misame_unscaled=p_misame_unscaled,
       p_misame_unscaled=p_misame_unscaled,
       p_vital_unscaled=p_vital_unscaled),
        file= paste0(here::here(),"/figure-data/figure_s4_trajectory_plots.RDS")
)
