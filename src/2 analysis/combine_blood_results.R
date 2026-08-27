# =============================================================================
# combine_blood_results.R
#
# Combines the separate blood-compartment clean result files into ONE tidy master
# table, analogous to the main milk results
# (adjusted_combined_arms_intervention_effects_results_clean.RDS).
#
# Binds every available blood clean file across {arm coding} x {adjustment}, tags
# each row, and harmonizes the schema toward the milk results (adds study, omic,
# tissue, label_f). Picks up adjusted files automatically once they exist.
#
# Output:
#   results/blood_compartment_all_intervention_effects_results_clean.RDS  (full)
#   results/blood_compartment_all_FDRsig_ATE.csv                          (slim: FDR-sig ATEs)
# =============================================================================

suppressMessages({library(dplyr)})
r <- paste0(here::here(), "/results/")

# Every blood clean file: tag = "<adjust><coding>"; path built from the two tags.
specs <- list(
  list(adj = "unadjusted", coding = "stratified", file = "blood_compartment_intervention_effects_results_clean.RDS"),
  list(adj = "unadjusted", coding = "combined",   file = "blood_compartment_combined_arms_intervention_effects_results_clean.RDS"),
  list(adj = "adjusted",   coding = "stratified", file = "blood_compartment_adjusted_intervention_effects_results_clean.RDS"),
  list(adj = "adjusted",   coding = "combined",   file = "blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")
)

read_one <- function(s) {
  p <- paste0(r, s$file)
  if (!file.exists(p)) { message("skip (not found): ", s$file); return(NULL) }
  readRDS(p) %>% mutate(adjustment = s$adj, arm_coding = s$coding)
}

blood <- bind_rows(lapply(specs, read_one))
if (nrow(blood) == 0) stop("No blood clean files found.")

# Harmonize toward the milk schema: study, omic, tissue, label_f.
blood <- blood %>%
  mutate(
    study = "MISAME-3",
    omic  = ifelse(grepl("Proteomics", dataset), "proteomics", "metabolomics"),
    tissue = dplyr::recode(dataset,
      MaternalPlasma        = "maternal plasma",
      VamsPrenatal          = "maternal VAMS (prenatal)",
      VamsPostnatalMaternal = "maternal VAMS (postnatal)",
      VamsPostnatalInfant   = "infant VAMS (postnatal)",
      ProteomicsDepleted    = "maternal blood (depleted)",
      ProteomicsNaive       = "maternal blood (naive)",
      .default = dataset),
    label_f = biomarker) %>%
  select(study, omic, tissue, dataset, compartment, visit, arm_coding, adjustment,
         contrast, measure, biomarker, label_f, est, cil, ciu,
         pval, pval_adj, pval_adj_global, sig, sigFDR, studytime)

saveRDS(blood, paste0(r, "blood_compartment_all_intervention_effects_results_clean.RDS"))

# Slim, portable CSV: FDR-significant ATEs only.
slim <- blood %>% filter(measure == "ATE", sigFDR == 1) %>%
  arrange(adjustment, arm_coding, dataset, visit, pval_adj)
write.csv(slim, paste0(r, "blood_compartment_all_FDRsig_ATE.csv"), row.names = FALSE)

cat("Combined blood results:\n")
cat("  total rows:", nrow(blood), "\n")
print(blood %>% count(adjustment, arm_coding, omic))
cat("\nFDR-significant ATEs (slim CSV):", nrow(slim), "rows across",
    dplyr::n_distinct(slim$dataset), "datasets\n")
