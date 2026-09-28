# =============================================================================
# 12-blood-compartment-intervention-effects.R   (arm-stratified)
#
# Estimates per-feature BEP average treatment effects (bioTMLE) in the MISAME-III
# blood compartments, separately by visit, for the Methods' cross-compartment
# blood analyses. Arms are the full 4-level factorial (Control, IFA-BEP, BEP-IFA,
# BEP-BEP), each contrasted with Control, mirroring the arm-stratified milk script
# 2_adjusted_analysis.R. Compartments: maternal plasma, prenatal VAMS, postnatal
# VAMS split by dyad into maternal and infant, and maternal blood proteomics
# (depleted and naive). FDR is applied later, in clean_blood_results.R. The
# covariate-adjusted stratified results feed script 20.
#
# Inputs : data/blood/merged_blood_datasets.RDS (1 data prep/7-blood-compartment-prep.R),
#          metadata/blood_component.Rdata
# Output : results/blood_compartment_[adjusted_]intervention_effects_results.RDS
# Settings (set before sourcing): BLOOD_ADJUST, BLOOD_WORKERS, BLOOD_SL_LIB_METAB,
#          BLOOD_SL_LIB_PROT, N_FEATURES_SUBSET (NULL = all features).
# [needs restricted data]
# =============================================================================

# Standalone runs load the config here; the orchestrators load it once and set BLOOD_ORCHESTRATED.
if (!exists("BLOOD_ORCHESTRATED")) { rm(list = ls()); source(paste0(here::here(), "/src/0-config.R")) }

load(file = paste0(here::here(), "/metadata/blood_component.Rdata"))
blood <- readRDS(paste0(here::here(), "/data/blood/merged_blood_datasets.RDS"))

if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- FALSE
# Cap parallel workers: SnowParam copies the data to each worker (Windows PSOCK),
# so an unbounded worker count runs out of memory on the 38k-feature VAMS.
if (!exists("BLOOD_WORKERS")) BLOOD_WORKERS <- 4
# SuperLearner library per compartment. The full milk library
# c("SL.mean","SL.glm","SL.glmnet","SL.xgboost") is about 31x too slow for the
# 19k-38k untargeted metabolite features, so both default to SL.glm; set
# BLOOD_SL_LIB_PROT before sourcing to use a larger library on proteomics.
if (!exists("BLOOD_SL_LIB_METAB")) BLOOD_SL_LIB_METAB <- c("SL.glm")
if (!exists("BLOOD_SL_LIB_PROT"))  BLOOD_SL_LIB_PROT  <- c("SL.glm")
sl_for <- function(label) if (grepl("Proteomics", label)) BLOOD_SL_LIB_PROT else BLOOD_SL_LIB_METAB
SCALE_IN_TMLE     <- TRUE               # standardize features inside bioTMLE (effects in SD units)
Wvars_blood       <- if (BLOOD_ADJUST) Wvars else c("arm", "dummy")  # adjusted = milk Wvars
.tag              <- if (BLOOD_ADJUST) "adjusted_" else ""            # output filename tag
if (!exists("N_FEATURES_SUBSET")) N_FEATURES_SUBSET <- NULL  # all features by default; set a number (e.g. 100) before sourcing for a quick test run

# Compartment spec: label -> (data.frame, feature set). Postnatal VAMS is split by dyad.
vp <- blood$vams_postnatal
specs <- list(
  MaternalPlasma         = list(df = blood$maternal_plasma,             Y = blood_components$maternal_plasma),
  VamsPrenatal           = list(df = blood$vams_prenatal,               Y = blood_components$vams_prenatal),
  VamsPostnatalMaternal  = list(df = dplyr::filter(vp, dyad=="mother"), Y = blood_components$vams_postnatal),
  VamsPostnatalInfant    = list(df = dplyr::filter(vp, dyad=="infant"), Y = blood_components$vams_postnatal),
  ProteomicsDepleted     = list(df = blood$proteomics_depleted,         Y = blood_components$proteomics_depleted),
  ProteomicsNaive        = list(df = blood$proteomics_naive,            Y = blood_components$proteomics_naive))

run_blood_compartment <- function(df, Yvars, label) {
  df <- df %>% filter(!is.na(arm), !is.na(visit))
  if (nrow(df) == 0) { warning(sprintf("[%s] no usable rows -- skipped.", label)); return(NULL) }
  Yvars <- intersect(Yvars, colnames(df))
  if (!is.null(N_FEATURES_SUBSET)) Yvars <- head(Yvars, N_FEATURES_SUBSET)
  sl <- sl_for(label)
  res <- df %>%
    group_by(visit) %>%
    do(res = try(run_bioTMLE(
      d = ., Wvars = Wvars_blood, bppar.debug = TRUE,
      g_lib = sl, Q_lib = sl,
      Yvars = Yvars, scale = SCALE_IN_TMLE,
      bppar.type = BiocParallel::SnowParam(workers = BLOOD_WORKERS))))
  names(res$res) <- paste0(label, "-", res$visit)
  res
}

results <- lapply(names(specs), function(lab)
  run_blood_compartment(specs[[lab]]$df, specs[[lab]]$Y, lab))
names(results) <- names(specs)

saveRDS(results, file = paste0(here::here(),
  "/results/blood_compartment_", .tag, "intervention_effects_results.RDS"))
