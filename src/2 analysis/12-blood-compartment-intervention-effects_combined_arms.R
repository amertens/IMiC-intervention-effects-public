
# =============================================================================
# 12-blood-compartment-intervention-effects_combined_arms.R   [COMBINED ARMS]
#
# Combined-arms companion, mirroring 2_adjusted_analysis_combined_arms.R.
# Binary Control vs BEP from the prep's `class` column = exposure-period-aware
# collapse (prenatal samples -> prenatal BEP; postnatal -> postpartum BEP).
#
# Analyze BY VISIT (= timePoint); FDR per visit x dataset in clean_blood_results.R.
# Postnatal VAMS split by dyad -> maternal vs infant circulation.
# Output: results/blood_compartment_combined_arms_intervention_effects_results.RDS
# =============================================================================

if (!exists("BLOOD_ORCHESTRATED")) { rm(list = ls()); source(paste0(here::here(), "/src/0-config.R")) }

load(file = paste0(here::here(), "/metadata/blood_component.Rdata"))
blood <- readRDS(paste0(here::here(), "/data/blood/merged_blood_datasets.RDS"))

if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- FALSE
if (!exists("BLOOD_WORKERS")) BLOOD_WORKERS <- 4   # cap PSOCK workers; unbounded OOMs on 38k-feat VAMS
# Per-compartment SL library: full library for small proteomics, SL.glm for the
# 19k-38k untargeted metabolomics (full library ~31x too slow there).
if (!exists("BLOOD_SL_LIB_METAB")) BLOOD_SL_LIB_METAB <- c("SL.glm")
if (!exists("BLOOD_SL_LIB_PROT"))  BLOOD_SL_LIB_PROT  <- c("SL.glm")
sl_for <- function(label) if (grepl("Proteomics", label)) BLOOD_SL_LIB_PROT else BLOOD_SL_LIB_METAB
SCALE_IN_TMLE     <- TRUE
Wvars_blood       <- if (BLOOD_ADJUST) Wvars else c("arm", "dummy")
.tag              <- if (BLOOD_ADJUST) "adjusted_" else ""
if (!exists("N_FEATURES_SUBSET")) N_FEATURES_SUBSET <- 100  # setup: first N feats/dataset; NULL=all

vp <- blood$vams_postnatal
specs <- list(
  MaternalPlasma         = list(df = blood$maternal_plasma,             Y = blood_components$maternal_plasma),
  VamsPrenatal           = list(df = blood$vams_prenatal,               Y = blood_components$vams_prenatal),
  VamsPostnatalMaternal  = list(df = dplyr::filter(vp, dyad=="mother"), Y = blood_components$vams_postnatal),
  VamsPostnatalInfant    = list(df = dplyr::filter(vp, dyad=="infant"), Y = blood_components$vams_postnatal),
  ProteomicsDepleted     = list(df = blood$proteomics_depleted,         Y = blood_components$proteomics_depleted),
  ProteomicsNaive        = list(df = blood$proteomics_naive,            Y = blood_components$proteomics_naive))

run_blood_compartment <- function(df, Yvars, label) {
  df <- df %>%
    mutate(arm = factor(ifelse(class == 1, "BEP", "Control"), levels = c("Control", "BEP"))) %>%
    filter(!is.na(arm), !is.na(visit))
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
  "/results/blood_compartment_", .tag, "combined_arms_intervention_effects_results.RDS"))
