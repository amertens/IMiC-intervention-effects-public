
# =============================================================================
# 7-blood-compartment-prep.R
#
# Builds blood-side analysis datasets for the Science-revision cross-compartment
# ATE, following Lishi's authoritative prep (metabolomProteomicsDataPath.Rmd).
# Linkage validated on real data (idBiospe match 93-100% per compartment).
#
# Compartments built:
#   maternal_plasma     prenatalPlasmaImputedCappedScaled.csv     (samples x feat)
#   vams_prenatal       prenatalVamsImputedCappedScaled.csv       (samples x feat, maternal)
#   vams_postnatal      postnatalVamsImputedCappedScaled.csv      (feat x samples -> t(); mother+infant)
#   proteomics_depleted ProteomicsDepletedMSStatsAndRFImputed.csv (samples x feat)  [primary]
#   proteomics_naive    proteomicsNaiveMSStatsAndRFImputed.csv    (samples x feat)  [supplementary]
#
# Treatment crosswalk: parsed idBiospe -> metadata_for_sharing.dta (idbs). The
# `code_bep_n` label encodes the full 2x2 factorial ("pre:..; post:.."), parsed
# into bep_prenatal / bep_postpartum.
#
# Output:
#   data/blood/merged_blood_datasets.RDS
#   metadata/blood_component.Rdata   (blood_components)
# =============================================================================

if (!exists("BLOOD_ORCHESTRATED")) { rm(list = ls()); source(paste0(here::here(), "/src/0-config.R")) }
library(haven)

dir_add  <- paste0(here::here(), "/data/additional datasets/")
meta_dta <- paste0(dir_add, "metadata_for_sharing.dta")

# ----------------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------------
# Feed the Imputed/Capped/Scaled prep files. run_bioTMLE(scale = TRUE) re-standardizes
# within each per-visit analysis subset, so ATEs are in WITHIN-VISIT SD units, not global
# SD (re-scaling a globally-scaled vector on a subset is NOT idempotent). This is fine for
# direction-only cross-compartment comparisons; do not compare SD magnitudes across
# datasets. scale=TRUE is kept to match the milk pipeline convention.
#
# Adjustment: UNADJUSTED (Wvars = c("arm","dummy")) by default. Set BLOOD_ADJUST=TRUE
# (the run_blood_full_adjusted.R orchestrator does this) to attach the milk-pipeline
# Wvars via the milk-ID bridge: blood idBiospe == milk BMID number -> subjid (1:1,
# verified) -> clean_baseline_covariates.RDS. (Override-aware so orchestrators preset it.)
if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- FALSE

# Feature subset: NULL = ALL features (default; correct, though the ~38k-feature VAMS
# takes hours). Set SETUP_N_FEATURES to a number (e.g. 100) BEFORE sourcing for a fast
# setup/test pass. Default is NULL so a standalone run never silently truncates features.
if (!exists("SETUP_N_FEATURES")) SETUP_N_FEATURES <- NULL

# In-utero exposure visits (class = prenatal BEP): enrollment, trimester 3, and BIRTH.
# Infant blood at delivery (acco) reflects PRENATAL exposure, not postpartum BEP (which
# has not started), so acco must key `class` on bep_prenatal, not bep_postpartum.
PRENATAL_TOKENS <- c("incl", "tri3", "acco")

# ----------------------------------------------------------------------------
# Treatment crosswalk: idBiospe -> factorial arm
# ----------------------------------------------------------------------------
if (!file.exists(meta_dta)) stop("Missing crosswalk: ", meta_dta)

arm_key <- read_dta(meta_dta) %>%
  mutate(idBiospe  = as.numeric(idbs),
         code_chr  = as.character(haven::as_factor(code_bep_n))) %>%
  filter(!is.na(idBiospe)) %>%
  distinct(idBiospe, code_chr) %>%
  mutate(
    bep_prenatal   = as.integer(grepl("pre:BEP",  code_chr, fixed = TRUE)),
    bep_postpartum = as.integer(grepl("post:BEP", code_chr, fixed = TRUE)),
    arm = factor(dplyr::case_when(
      bep_prenatal == 0 & bep_postpartum == 0 ~ "Control",
      bep_prenatal == 0 & bep_postpartum == 1 ~ "IFA/BEP",
      bep_prenatal == 1 & bep_postpartum == 0 ~ "BEP/IFA",
      bep_prenatal == 1 & bep_postpartum == 1 ~ "BEP/BEP"),
      levels = c("Control", "IFA/BEP", "BEP/IFA", "BEP/BEP")))

#' Period-aware binary treatment (= combined-arms `class`): prenatal samples keyed
#' on prenatal BEP, postnatal on postpartum BEP. Robust to files that mix periods.
add_class_binary <- function(df) {
  prenatal <- df$timePoint %in% PRENATAL_TOKENS
  df$class <- ifelse(prenatal, df$bep_prenatal, df$bep_postpartum)
  df
}

# Adjusted-covariate inputs (built once, only when adjusting).
# Milk-ID bridge: milk bmid = "{timepoint}_{number}", and that number == blood
# idBiospe (verified 290/290, 1:1 to subjid). So idBiospe -> subjid -> Wvars.
if (BLOOD_ADJUST) {
  .milk <- readRDS(paste0(here::here(), "/data/merged_analysis_datasets.RDS"))
  idbiospe_to_subjid <- .milk %>%
    filter(study == "Misame") %>%
    transmute(idBiospe = as.numeric(sub(".*[_-]", "", bmid)),
              subjid    = as.integer(subjid)) %>%
    distinct() %>% filter(!is.na(idBiospe), !is.na(subjid))
  .cov_vars <- setdiff(Wvars, "arm")   # bring all milk Wvars EXCEPT treatment arm
  baseline_cov <- readRDS(paste0(here::here(), "/data/clean_baseline_covariates.RDS")) %>%
    filter(studyid == "MISAME-3") %>%
    mutate(subjid = as.integer(subjid)) %>%
    select(subjid, all_of(.cov_vars)) %>%
    distinct(subjid, .keep_all = TRUE)
}

attach_covariates <- function(df) {
  df$dummy <- 1L                       # always present so unadjusted models still run
  if (!BLOOD_ADJUST) return(df)
  # idBiospe -> subjid -> Wvars. NOTE: adjusted models are RESTRICTED to milk-linked dyads
  # (samples whose mother has a milk sample); samples without a milk/baseline match have no
  # covariate link and are dropped here. This is a selected subset (milk availability is a
  # post-randomization characteristic) -- log how many are dropped so the restriction is
  # visible, and report it alongside the adjusted results.
  out <- df %>%
    left_join(idbiospe_to_subjid, by = "idBiospe") %>%
    left_join(baseline_cov, by = "subjid")
  n_drop <- sum(is.na(out$subjid))
  if (n_drop > 0) message(sprintf("  [attach_covariates] adjusted: dropped %d/%d samples with no milk-linked covariates",
                                  n_drop, nrow(out)))
  out %>% filter(!is.na(subjid))
}

attach_treatment <- function(df, label) {
  out <- df %>%
    left_join(arm_key, by = "idBiospe") %>%
    add_class_binary() %>%
    attach_covariates()
  message(sprintf("[%s] n=%d  treatment match=%.1f%%  visits={%s}",
                  label, nrow(out), 100 * mean(!is.na(out$arm)),
                  paste(sort(unique(out$timePoint)), collapse = ",")))
  out %>% mutate(visit = timePoint) %>% filter(!is.na(class))
}

# ----------------------------------------------------------------------------
# Sample-id parsers (strsplit-based; return timePoint, idBiospe [, dyad])
# ----------------------------------------------------------------------------
nth <- function(lst, i) vapply(lst, function(z) if (length(z) >= i) z[i] else NA_character_, "")

parse_plasma <- function(x) {                 # "plasma;tri3_sapl_618;05/03/2021"
  p <- strsplit(nth(strsplit(x, ";"), 2), "_")
  tibble(timePoint = nth(p, 1), idBiospe = as.numeric(nth(p, 3)))
}
parse_vams_prenatal <- function(x) {          # "VAMS;incl_sa10ab_512;711732"
  p <- strsplit(nth(strsplit(x, ";"), 2), "_")
  tibble(timePoint = nth(p, 1), idBiospe = as.numeric(nth(p, 3)), dyad = "mother")
}
parse_vams_postnatal <- function(x) {         # "pn12_sa10ab_302_e_711604" | "..._302_711604"
  p <- strsplit(x, "_")
  tibble(timePoint = nth(p, 1), idBiospe = as.numeric(nth(p, 3)),
         dyad = ifelse(nth(p, 4) == "e", "infant", "mother"))
}
parse_proteomics <- function(x) {             # "pn12_106"
  p <- strsplit(gsub('"', "", x), "_")
  tibble(timePoint = nth(p, 1), idBiospe = as.numeric(nth(p, 2)))
}

# ----------------------------------------------------------------------------
# Compartment loaders (data.table::fread -- fast on wide files)
# ----------------------------------------------------------------------------
# samples x features, sample id in column 1
load_rows <- function(file, parser, label) {
  sel <- if (!is.null(SETUP_N_FEATURES)) seq_len(SETUP_N_FEATURES + 1) else NULL
  m <- data.table::fread(paste0(dir_add, file), data.table = FALSE, showProgress = FALSE, select = sel)
  feat <- m[, -1, drop = FALSE]
  list(df = bind_cols(parser(m[[1]]), feat) %>% attach_treatment(label), feats = colnames(feat))
}
# features x samples, sample id in header -> transpose
load_cols <- function(file, parser, label) {
  nr <- if (!is.null(SETUP_N_FEATURES)) SETUP_N_FEATURES else Inf
  m <- data.table::fread(paste0(dir_add, file), data.table = FALSE, showProgress = FALSE, nrows = nr)
  featnames <- m[[1]]; samples <- colnames(m)[-1]
  mat <- t(as.matrix(m[, -1, drop = FALSE])); rownames(mat) <- samples; colnames(mat) <- featnames
  list(df = bind_cols(parser(samples), as.data.frame(mat, check.names = FALSE)) %>% attach_treatment(label),
       feats = featnames)
}

maternal_plasma <- load_rows("prenatalPlasmaImputedCappedScaled.csv", parse_plasma,        "maternal_plasma")
vams_prenatal   <- load_rows("prenatalVamsImputedCappedScaled.csv",   parse_vams_prenatal, "vams_prenatal")
vams_postnatal  <- load_cols("postnatalVamsImputedCappedScaled.csv",  parse_vams_postnatal,"vams_postnatal")
proteomics_dep  <- load_rows("ProteomicsDepletedMSStatsAndRFImputed.csv", parse_proteomics, "proteomics_depleted")
proteomics_nai  <- load_rows("proteomicsNaiveMSStatsAndRFImputed.csv",    parse_proteomics, "proteomics_naive")

# ----------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------
blood_components <- list(
  maternal_plasma     = maternal_plasma$feats,
  vams_prenatal       = vams_prenatal$feats,
  vams_postnatal      = vams_postnatal$feats,
  proteomics_depleted = proteomics_dep$feats,
  proteomics_naive    = proteomics_nai$feats)

dir.create(paste0(here::here(), "/data/blood"), showWarnings = FALSE)
saveRDS(list(
  maternal_plasma     = maternal_plasma$df,
  vams_prenatal       = vams_prenatal$df,
  vams_postnatal      = vams_postnatal$df,   # split by dyad downstream (mother/infant)
  proteomics_depleted = proteomics_dep$df,
  proteomics_naive    = proteomics_nai$df),
  file = paste0(here::here(), "/data/blood/merged_blood_datasets.RDS"))
save(blood_components, file = paste0(here::here(), "/metadata/blood_component.Rdata"))

message("Done. arm = 4-level factorial (stratified); class = period-aware binary (combined).")
