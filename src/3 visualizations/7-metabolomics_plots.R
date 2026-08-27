# =============================================================================
# src/3 visualizations/7-metabolomics_plots.R
#
# Reads:  metadata/milk_component.Rdata
#         results/adjusted_intervention_effects_res_untargeted_metabolomics_clean.RDS
# Writes: figures/volcano_plot_metabolomics.png
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
res <- readRDS(paste0(here::here(),"/results/adjusted_intervention_effects_res_untargeted_metabolomics_clean.RDS"))
head(res)

# res <- Map(extract_res, d$res) %>% rbindlist(., idcol='studytime') %>% 
#   mutate(outcome_group='untargeted metabolomics', label_f=biomarker ) %>% 
#   as.data.frame()
# 
# res <- res %>% mutate(
#   studytime=as.character(studytime),
#   studytime = case_when(
#     studytime=="Misame-1" ~ "Misame (14-21 days)",
#     studytime=="Misame-2" ~ "Misame (1-2 mo.)",
#     studytime=="Misame-3" ~ "Misame (3-4 mo.)",
#     studytime=="Vital-40" ~ "Vital (1.5 mo.)",
#     studytime=="Vital-56" ~ "Vital (2 mo.)",
#     studytime=="Elicit-1" ~ "Elicit (1 mo.)",
#     studytime=="Elicit-5" ~ "Elicit (5 mo.)"
#   ),
#   studytime=factor(studytime, levels=c("Misame (14-21 days)", "Misame (1-2 mo.)", "Misame (3-4 mo.)",
#                                        "Vital (1.5 mo.)", "Vital (2 mo.)", "Elicit (1 mo.)", "Elicit (5 mo.)"))
# )
# 
# Earlier in the script we already loaded the cleaned `res` data.frame from
# `adjusted_intervention_effects_res_untargeted_metabolomics_clean.RDS`, so we
# don't need to re-extract from the raw list form (which lived in the
# misspelled-and-now-removed `metabalomics_intervention_effects_results.RDS`).
# Skip the re-extraction; the loaded `res` already has all required columns.

# (legacy block left for reference)
# d <- readRDS(paste0(here::here(),"/results/metabolomics_intervention_effects_results.RDS"))
# res <- Map(extract_res, d$res) %>% rbindlist(., idcol='studytime') %>%
#   mutate(outcome_group='untargeted metabolomics', label_f=biomarker ) %>%
#   as.data.frame()

res <- res %>% mutate(
  studytime=as.character(studytime),
  studytime = case_when(
    studytime=="Misame-1" ~ "Misame (14-21 days)",
    studytime=="Misame-2" ~ "Misame (1-2 mo.)",
    studytime=="Misame-3" ~ "Misame (3-4 mo.)",
    studytime=="Vital-40" ~ "Vital (1.5 mo.)",
    studytime=="Vital-56" ~ "Vital (2 mo.)",
    studytime=="Elicit-1" ~ "Elicit (1 mo.)",
    studytime=="Elicit-5" ~ "Elicit (5 mo.)",
    # The cleaned results already store studytime in readable form
    # ("Misame (14-21 days)", etc.); pass those through unchanged. Without this
    # fallback, case_when() returns NA for already-formatted values, which empties
    # every per-study volcano panel.
    TRUE ~ studytime
  ),
  studytime=factor(studytime, levels=c("Misame (14-21 days)", "Misame (1-2 mo.)", "Misame (3-4 mo.)",
                                       "Vital (1.5 mo.)", "Vital (2 mo.)", "Elicit (1 mo.)", "Elicit (5 mo.)"))
)


res <- res %>% filter(measure=="ATE", biomarker!="ARM") %>%
  group_by(studytime, biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup() %>% arrange(chi_pval_adj)



p_metabolomics_plot1 <- plot_imic_volcano(res %>% filter(grepl("Elicit",studytime)), labels=10)
p_metabolomics_plot2 <- plot_imic_volcano(res %>% filter(grepl("Vital",studytime)),  labels=10)
p_metabolomics_plot3 <- plot_imic_volcano(res %>% filter(grepl("Misame",studytime)), labels=10)

# Defensive: combine + save in tryCatch so an empty per-study filter (which
# returns a stub plot from plot_imic_volcano) doesn't crash cowplot's gtable
# alignment for the combined grid.
tryCatch({
  p_metabolomics_plot <- cowplot::plot_grid(
    p_metabolomics_plot3, p_metabolomics_plot2, p_metabolomics_plot1,
    ncol = 1, rel_heights = c(1, 1, 1))
  ggsave(p_metabolomics_plot,
         file = paste0(here::here(), "/figures/volcano_plot_metabolomics.png"),
         width = 18, height = 10)
}, error = function(e) {
  message("WARNING: skipped combined metabolomics volcano (", conditionMessage(e), ")")
  # fall back to saving individual panels
  for (i in seq_len(3)) {
    plt <- list(p_metabolomics_plot3, p_metabolomics_plot2, p_metabolomics_plot1)[[i]]
    nm  <- c("misame","vital","elicit")[i]
    tryCatch(
      ggsave(plt, file = paste0(here::here(),"/figures/volcano_plot_metabolomics_",nm,".png"),
             width = 18, height = 4),
      error = function(e2) message("  also skipped ", nm, ": ", conditionMessage(e2))
    )
  }
})


#to do: get the metadata labels for the metabolites. Look for commonalities across the studies


