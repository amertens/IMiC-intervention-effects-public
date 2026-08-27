# =============================================================================
# src/3 visualizations/triglyceride_trajectories.R
#
# Reads:  data/merged_analysis_datasets.RDS
#         metadata/milk_component.Rdata
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
dfull <- readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

d_metabolomics <- dfull %>% filter(study=="Misame") %>%
  select("visit","arm","mbmi", misame_milk_components$metabolomics[grepl("tg.",misame_milk_components$metabolomics)])

d_metabolomics <- d_metabolomics %>% gather(key="metabolite", value="value", -c(visit,arm,mbmi))

#drop metabolites with 2 values (rare ones coded as abundance)
d_metabolomics <- d_metabolomics %>% group_by(metabolite) %>% 
  mutate(n_distinct=n_distinct(value), value=scale(value)) %>% ungroup() %>%
  filter(n_distinct>2)

#faceted line ggplot with metabolites over time
p_metabolomics_control <- ggplot(d_metabolomics %>% filter(arm=="Control"), aes(x=visit, y=value, group=metabolite), alpha=0.25) + 
  stat_smooth(geom='line', alpha=0.1, se=FALSE)+theme_minimal()
p_metabolomics_control


#faceted line ggplot with metabolites over time
p_metabolomics<- ggplot(d_metabolomics, aes(x=visit, y=value, group=metabolite), alpha=0.25) + 
  facet_wrap(~arm, scales="fixed", nrow=1) +
  stat_smooth(geom='line', alpha=0.1, se=FALSE)+theme_minimal()
p_metabolomics

p_metabolomics_strat <- ggplot(d_metabolomics, aes(x=visit, y=value, group=metabolite), alpha=0.25) + 
  facet_wrap(mbmi < 18.5~arm, scales="fixed", nrow=2) +
  stat_smooth(geom='line', alpha=0.1, se=FALSE)+theme_minimal()
p_metabolomics_strat

p_metabolomics <- ggplot(d_metabolomics, aes(x=visit, y=value, group=metabolite, color=mbmi), alpha=0.25) + 
  facet_wrap(~arm, scales="fixed", nrow=2) +
  stat_smooth(aes(color=mbmi),geom='line', alpha=0.1, se=FALSE)+theme_minimal()
p_metabolomics

#need to fix to show change in relationshop to BMI
p_metabolomics <- ggplot(d_metabolomics, 
                         aes(x=visit, y=value, group=metabolite, color=mbmi, fill=mbmi), alpha=0.25) + 
  facet_wrap(~arm, scales="fixed", nrow=2) +
  geom_point() + geom_line() + theme(legend.position="bottom")
p_metabolomics

d_metabolomics <- dfull %>% filter(study=="Vital", arm=="Control") %>%
  select("visit", misame_milk_components$metabolomics[grepl("tg.",misame_milk_components$metabolomics)])

d_metabolomics <- d_metabolomics %>% gather(key="metabolite", value="value", -visit)

#drop metabolites with 2 values (rare ones coded as abundance)
d_metabolomics <- d_metabolomics %>% group_by(metabolite) %>% 
  mutate(n_distinct=n_distinct(value), value=scale(value)) %>% ungroup() %>%
  filter(n_distinct>2)

p_vital <- ggplot(d_metabolomics, aes(x=visit, y=value, group=metabolite), alpha=0.25) + 
  geom_line(alpha=0.1)+theme_minimal()
p_vital


