#-------------------------------------------------------------------------------
# Build Table S7: native-unit means (± SD) + ATEs for primary + secondary outcomes
# per intervention arm × visit × trial.
#
# Output: results/tables/table_s7_primary_secondary_native_units.csv (long)
#         results/tables/table_s7_primary_secondary_native_units_wide.csv (wide)
#
# Reviewer 2 §2.4 requested absolute concentrations alongside the z-scored ATEs in
# Fig. 3. This script consumes the pre-existing per-arm-and-visit unscaled effect
# estimates and joins them with N + SD computed from the raw merged dataset.
#
# Owner: Andrew. Effort: ~10 min once the source files are confirmed.
#-------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

RESULTS  <- paste0(here::here(), "/results/")
DATA     <- paste0(here::here(), "/data/")
OUT      <- paste0(here::here(), "/results/tables/")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

#-- 1. Load inputs ----------------------------------------------------------
unscaled <- readRDS(file.path(RESULTS, "adjusted_combined_arms_intervention_effects_unscaled_results_clean.RDS"))
raw      <- readRDS(file.path(DATA,    "merged_analysis_datasets.RDS"))

cat("Loaded unscaled results: ", nrow(unscaled), " rows\n", sep="")
cat("Loaded raw analytic data: ", nrow(raw), " rows, ", ncol(raw), " cols\n\n", sep="")

#-- 2. Define primary + secondary outcome categories ------------------------
# Primary  = macronutrients + micronutrients + B-vitamins + fat-soluble vitamins
# Secondary = HMOs + bioactive proteins
primary_cats <- c(
  "Macronutrient",
  "Micronutrient",
  "B1", "B2", "B3 or related", "B6", "Other B vitamins"
)
secondary_cats <- c(
  "Bioactive",
  "Individual HMO containing Fucose",
  "Individual HMO containing Sialic Acid",
  "Other individual HMO"
)
keep_cats <- c(primary_cats, secondary_cats)

# Vitamin A and the tocopherols come through `unscaled` with category/label == NA
# (a biomarker-name-casing miss against the upstream category lookup, same failure
# mode as the known "vitamin.a" join gotcha), which silently dropped them from every
# prior version of this table even though they are primary micronutrient outcomes
# reported in the main text. Patch category/label for these three by name before
# filtering rather than dropping them.
# Fgf.21/Iga come through as category "Other individual HMO" (an upstream placeholder-
# category recode), contradicting the manuscript's own Methods text, which names both as
# example "selected bioactive proteins" alongside secretory IgA, calprotectin, leptin,
# insulin. Fsh/Lh/Diversity/Evenness come through tagged "Bioactive" even though FSH/LH
# are never described as an outcome anywhere in the manuscript and Diversity/Evenness are
# the microbiome alpha-diversity indices the manuscript explicitly calls "exploratory," not
# secondary. Re-categorize all six to match the manuscript text (see
# functions/data_cleaning_functions.R for the same fix at the source; patched again here
# because this RDS was generated before that fix and a full pipeline rebuild is out of
# scope for this table alone).
unscaled <- unscaled %>%
  mutate(
    category = case_when(
      biomarker %in% c("Vitamin.a", "A.tocopherol", "G.tocopherol") ~ "Micronutrient",
      biomarker %in% c("Fgf.21", "Iga") ~ "Bioactive",
      biomarker %in% c("Fsh", "Lh") ~ NA_character_,
      biomarker %in% c("Diversity", "Evenness") ~ "Exploratory (microbiome)",
      TRUE ~ category
    ),
    label = case_when(
      biomarker == "Vitamin.a"    ~ "Vitamin A",
      biomarker == "A.tocopherol" ~ "α-Tocopherol",
      biomarker == "G.tocopherol" ~ "γ-Tocopherol",
      biomarker == "Fgf.21"       ~ "FGF-21",
      biomarker == "Iga"          ~ "Secretory IgA",
      TRUE ~ label
    )
  )

dat <- unscaled %>%
  ungroup() %>%
  filter(category %in% keep_cats) %>%
  mutate(
    outcome_class = ifelse(category %in% primary_cats,   "primary",
                    ifelse(category %in% secondary_cats, "secondary", NA))
  )
cat("Filtered to primary + secondary outcomes: ", nrow(dat), " rows\n", sep="")
cat("Categories retained:\n"); print(table(dat$category, dat$outcome_class))

#-- 3. Split per-arm means from ATE rows ------------------------------------
# `measure == "MN"` rows give per-arm means (cil/ciu = 95% CI for the mean)
# `measure == "ATE"` rows give intervention effects (cil/ciu = 95% CI for the ATE)
means <- dat %>%
  filter(measure == "MN") %>%
  transmute(study, visit, biomarker, label, category, outcome_class,
            arm   = contrast,
            mean  = est,
            mean_cil = cil,
            mean_ciu = ciu)

ates  <- dat %>%
  filter(measure == "ATE") %>%
  transmute(study, visit, biomarker, label, category, outcome_class,
            contrast,
            ate    = est,
            ate_cil = cil,
            ate_ciu = ciu,
            pval, pval_adj, sigFDR)

#-- 4. Compute N and SD per arm × study × visit × biomarker from raw data ---
# The raw `arm` column has trial-specific arm labels (e.g., "Control", "BEP/BEP",
# "Nico", "BEP+ExBf", "BEP+ExBf+AZT", "Nico+Az.", "Az.", "IFA/BEP", "BEP/IFA").
# The combined-arms analysis groups multiple raw arms into a single "BEP" or
# "Nico" contrast level. Replicate that grouping here.
trial_visit_lookup <- function(study, visit_num) {
  # Map raw numeric visits → result-table visit labels
  # MISAME visits: 1, 2, 3 → "14--21 days", "1-2 mo.", "3-4 mo."
  # Mumta-LW (Vital) visits: 40, 56 → "1.5 mo.", "2 mo." (days postpartum approx)
  # ELICIT visits: 1, 5 → "1 mo.", "5 mo."
  if (study == "Misame") {
    return(c(`1` = "14-21 days", `2` = "1-2 mo.", `3` = "3-4 mo.")[as.character(visit_num)])
  } else if (study == "Vital") {
    return(c(`40` = "1.5 mo.", `56` = "2 mo.")[as.character(visit_num)])
  } else if (study == "Elicit") {
    return(c(`1` = "1 mo.", `5` = "5 mo.")[as.character(visit_num)])
  }
  NA_character_
}

# Build an arm-grouping rule consistent with the combined-arms intervention model
arm_to_contrast <- function(arm_raw, study) {
  arm_raw <- as.character(arm_raw)
  # Control levels
  if (arm_raw == "Control") return("Control")
  if (study == "Misame" && arm_raw == "IFA/BEP")     return("Control")  # no postnatal BEP
  if (study == "Misame" && arm_raw == "BEP/IFA")     return("BEP")      # postnatal BEP only
  if (study == "Misame" && arm_raw == "BEP/BEP")     return("BEP")      # postnatal BEP (combined)
  if (study == "Vital"  && arm_raw %in% c("BEP+ExBf", "BEP+ExBf+AZT")) return("BEP")
  if (study == "Elicit" && arm_raw == "Az.")         return("Control")  # azithromycin-only
  if (study == "Elicit" && arm_raw %in% c("Nico", "Nico+Az.")) return("Nico")
  NA_character_
}

# Apply mapping
raw2 <- raw %>%
  mutate(
    visit_label = mapply(trial_visit_lookup, study, visit),
    arm_grouped = mapply(arm_to_contrast, arm, study)
  ) %>%
  filter(!is.na(arm_grouped))

cat("\nRaw rows after arm grouping: ", nrow(raw2), "\n", sep="")
cat("Arm grouping summary:\n"); print(table(raw2$study, raw2$arm_grouped, useNA = "ifany"))

#-- 5. For each biomarker in primary/secondary, compute N, mean, SD ---------
biomarker_codes <- unique(c(means$biomarker, ates$biomarker))
cat("\nUnique biomarker codes to look up in raw data: ", length(biomarker_codes), "\n", sep="")

# Some biomarker codes may not exist as column names (typical name-case differences).
# Build a case-insensitive lookup, with an explicit alias table for biomarkers
# whose raw-data column name differs from the results-table biomarker code.
# (E.g., "Total carbohydrate" in MISAME-III is column `carbohydrate`; ELICIT/Mumta-LW
# use `cho` for the same outcome, which already case-insensitively matches biomarker
# code "Cho".)
biomarker_aliases <- c(
  "Total carbohydrate" = "carbohydrate"
)

raw_cols <- colnames(raw2)
lookup_col <- function(code) {
  # Explicit alias takes precedence
  if (code %in% names(biomarker_aliases)) {
    alias <- biomarker_aliases[[code]]
    ci <- which(tolower(raw_cols) == tolower(alias))
    if (length(ci) == 1) return(raw_cols[ci])
  }
  ci <- which(tolower(raw_cols) == tolower(code))
  if (length(ci) == 1) return(raw_cols[ci])
  NA_character_
}

raw_descriptives <- list()
for (bm in biomarker_codes) {
  col <- lookup_col(bm)
  if (is.na(col)) next
  sub <- raw2 %>%
    filter(!is.na(.data[[col]])) %>%
    group_by(study, visit_label, arm_grouped) %>%
    summarise(
      N   = n(),
      raw_mean  = mean(.data[[col]], na.rm = TRUE),
      raw_sd    = sd(.data[[col]],   na.rm = TRUE),
      raw_median = median(.data[[col]], na.rm = TRUE),
      raw_q25    = quantile(.data[[col]], 0.25, na.rm = TRUE),
      raw_q75    = quantile(.data[[col]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(biomarker = bm)
  raw_descriptives[[bm]] <- sub
}
raw_desc_df <- bind_rows(raw_descriptives)
cat("Biomarkers matched to raw data columns: ", length(raw_descriptives),
    " / ", length(biomarker_codes), "\n", sep="")

#-- 6. Build the long-form Table S7 -----------------------------------------
# Join means (per-arm) with raw descriptives (N, SD, median, IQR)
means_with_n <- means %>%
  left_join(raw_desc_df, by = c("study", "visit" = "visit_label", "biomarker", "arm" = "arm_grouped"))

# Ensure ATE row uses the intervention-arm label
ates_n <- ates %>%
  rename(arm = contrast) %>%
  mutate(measure = "ATE", mean = ate, mean_cil = ate_cil, mean_ciu = ate_ciu) %>%
  left_join(raw_desc_df %>% select(study, visit_label, biomarker, arm_grouped),
            by = c("study", "visit" = "visit_label", "biomarker", "arm" = "arm_grouped"))

#-- 7. Build the wide-form Table S7 (one row per biomarker × study × visit) -
wide <- means_with_n %>%
  select(outcome_class, category, label, biomarker, study, visit, arm,
         N, raw_mean, raw_sd, raw_median, raw_q25, raw_q75) %>%
  pivot_wider(
    names_from  = arm,
    values_from = c(N, raw_mean, raw_sd, raw_median, raw_q25, raw_q75),
    names_glue  = "{arm}_{.value}"
  ) %>%
  left_join(
    ates %>%
      transmute(study, visit, biomarker,
                ATE = ate, ATE_cil = ate_cil, ATE_ciu = ate_ciu,
                pval, pval_adj, sigFDR),
    by = c("study", "visit", "biomarker")
  ) %>%
  arrange(outcome_class, category, label, study, visit)

#-- 8. Write outputs --------------------------------------------------------
write.csv(means_with_n, file.path(OUT, "table_s7_primary_secondary_native_units.csv"),
          row.names = FALSE)
write.csv(wide, file.path(OUT, "table_s7_primary_secondary_native_units_wide.csv"),
          row.names = FALSE)

cat("\nWrote:\n  ", file.path(OUT, "table_s7_primary_secondary_native_units.csv"), "\n",
    "  ", file.path(OUT, "table_s7_primary_secondary_native_units_wide.csv"), "\n", sep="")
cat("\nLong rows: ", nrow(means_with_n), "\n",
    "Wide rows: ", nrow(wide), "\n", sep="")

#-- 9. Quick QC: how many biomarkers got matched? ---------------------------
cat("\nWide-form preview (head 8):\n")
print(head(wide, 8))

cat("\nSummary by outcome class:\n")
print(wide %>% count(outcome_class, category))

unmatched <- setdiff(biomarker_codes, names(raw_descriptives))
if (length(unmatched)) {
  cat("\nBiomarkers NOT found in raw data (will lack N/SD; effect estimates still present):\n")
  cat(paste0("  - ", head(unmatched, 30), collapse = "\n"), "\n", sep="")
  if (length(unmatched) > 30) cat("  ... (", length(unmatched), " total)\n", sep="")
}
