
# =============================================================================
# 12-blood-compartment-intervention-effects.R   [ARM-STRATIFIED]
#
# bioTMLE ATE on the blood compartments from 7-blood-compartment-prep.R,
# mirroring the arm-stratified milk script 2_adjusted_analysis.R.
#
# STRATIFIED = full 4-level MISAME-3 factorial arm (Control / IFA-BEP / BEP-IFA /
# BEP-BEP); bioTMLE contrasts each vs Control. Analyze BY VISIT (= timePoint);
# FDR per visit x dataset in clean_blood_results.R.
#
# Postnatal VAMS is split by dyad -> maternal vs infant circulation.
# Output: results/blood_compartment_intervention_effects_results.RDS
# =============================================================================

if (!exists("BLOOD_ORCHESTRATED")) { rm(list = ls()); source(paste0(here::here(), "/src/0-config.R")) }

load(file = paste0(here::here(), "/metadata/blood_component.Rdata"))
blood <- readRDS(paste0(here::here(), "/data/blood/merged_blood_datasets.RDS"))

if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- FALSE
# Cap parallel workers: SnowParam copies the data to EACH worker (Windows PSOCK),
# so unbounded workers OOM on the 38k-feature VAMS. Keep small. (override-aware)
if (!exists("BLOOD_WORKERS")) BLOOD_WORKERS <- 4
# Per-compartment SuperLearner library: the milk-analysis full library
# c("SL.mean","SL.glm","SL.glmnet","SL.xgboost") is feasible for the small
# proteomics datasets but ~31x too slow for 19k-38k untargeted metabolite features,
# so metabolomics defaults to SL.glm. Orchestrators override these.
if (!exists("BLOOD_SL_LIB_METAB")) BLOOD_SL_LIB_METAB <- c("SL.glm")
if (!exists("BLOOD_SL_LIB_PROT"))  BLOOD_SL_LIB_PROT  <- c("SL.glm")
sl_for <- function(label) if (grepl("Proteomics", label)) BLOOD_SL_LIB_PROT else BLOOD_SL_LIB_METAB
SCALE_IN_TMLE     <- TRUE               # per instruction
Wvars_blood       <- if (BLOOD_ADJUST) Wvars else c("arm", "dummy")  # adjusted = milk Wvars
.tag              <- if (BLOOD_ADJUST) "adjusted_" else ""            # output filename tag
if (!exists("N_FEATURES_SUBSET")) N_FEATURES_SUBSET <- NULL  # ALL features by default; set a number (e.g. 100) before sourcing for a fast setup pass

# Compartment spec: label -> (data.frame, feature set). Postnatal VAMS split by dyad.
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
