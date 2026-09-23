#-------------------------------------------------------------------------------
# Build the machine-readable companion of Table S1 for the TERTIARY outcomes:
# native-unit intervention effects for every targeted (MxP Quant 500) metabolite,
# including each triglyceride species, by trial x visit.
#
# Output: results/tables/targeted_metabolites_native_units.csv
#
# Why: round-2 Reviewer 1 (point 4e) asked where the native-unit triglyceride-species
# estimates are. Table S1 (formerly S7; build_table_s1.R) covers primary + secondary
# outcomes only, so this exports the same combined-arm, covariate-adjusted TMLE effects
# for the tertiary panel from the same upstream unscaled-effects run.
#
# One row per study x visit x metabolite x contrast. Columns: model-based marginal means
# for the control and intervention arms, the native-unit ATE with 95% CI, its raw P and
# BH q (native-unit model), and, for cross-reference with Fig. 5A / the text, the q-value
# from the standardized (Z-scored) model, which is the FDR the paper reports.
#-------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

RESULTS <- paste0(here::here(), "/results/")
OUT     <- paste0(here::here(), "/results/tables/")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

unscaled <- readRDS(file.path(RESULTS, "adjusted_combined_arms_intervention_effects_unscaled_results_clean.RDS"))
scaled   <- readRDS(file.path(RESULTS, "adjusted_combined_arms_intervention_effects_results_clean.RDS"))
stopifnot("tertiary" %in% unscaled$outcome_group)

study_lab <- c(Elicit = "ELICIT", Misame = "MISAME-III", Vital = "Mumta-LW")
u <- unscaled %>% filter(outcome_group == "tertiary")

means <- u %>% filter(measure == "MN") %>%
  mutate(arm = ifelse(contrast == "Control", "control_mean", "intervention_mean")) %>%
  select(study, visit, biomarker, arm, est) %>%
  distinct(study, visit, biomarker, arm, .keep_all = TRUE) %>%
  pivot_wider(names_from = arm, values_from = est)

ate <- u %>% filter(measure == "ATE") %>%
  select(study, visit, biomarker, label, category, contrast,
         ATE = est, ATE_cil = cil, ATE_ciu = ciu, pval_native = pval, qval_native = pval_adj)

q_z <- scaled %>% filter(measure == "ATE", outcome_group == "tertiary") %>%
  select(study, visit, biomarker, contrast, ATE_sd = est, qval_standardized = pval_adj)

out <- ate %>%
  left_join(means, by = c("study", "visit", "biomarker")) %>%
  left_join(q_z,   by = c("study", "visit", "biomarker", "contrast")) %>%
  mutate(study = unname(study_lab[study]),
         pct_change = 100 * ATE / control_mean) %>%
  select(study, visit, category, biomarker, label, contrast,
         control_mean, intervention_mean, ATE, ATE_cil, ATE_ciu, pct_change,
         pval_native, qval_native, ATE_sd, qval_standardized) %>%
  arrange(study, visit, category, biomarker)

f <- file.path(OUT, "targeted_metabolites_native_units.csv")
write.csv(out, f, row.names = FALSE)
cat("wrote", f, "|", nrow(out), "rows |", n_distinct(out$biomarker), "metabolites |",
    sum(out$category == "Triglycerides", na.rm = TRUE), "triglyceride rows |",
    sum(out$qval_standardized < 0.05, na.rm = TRUE), "FDR-significant (standardized model)\n")
