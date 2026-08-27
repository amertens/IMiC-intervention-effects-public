# =============================================================================
# 18-directional-mummichog.R
#
# Directional Mummichog: runs the UP and DOWN feature sets separately so each
# pathway's DIRECTION can be stated per compartment (mummichog enrichment is
# otherwise direction-agnostic). For each compartment x ion mode we run:
#   UP   : only features with est>0 are "significant" (p kept; others set to p=1)
#   DOWN : only features with est<0 are "significant"
# the full m/z list is the background in both cases.
#
# Compartments: Milk (adjusted-published, pooled visits), MaternalPlasma (pn12),
# MaternalVAMS (pn56), InfantVAMS (pooled postnatal). Adjusted blood.
#
# Output: results/mummichog_output_directional/, results/directional_pathways.csv
# =============================================================================

suppressMessages({library(dplyr); library(data.table)})
root <- paste0(here::here(), "/"); add <- paste0(root, "data/additional datasets/")
indir <- paste0(root, "results/mummichog_input_directional");  dir.create(indir, showWarnings = FALSE, recursive = TRUE)
outroot <- paste0(root, "results/mummichog_output_directional"); dir.create(outroot, showWarnings = FALSE, recursive = TRUE)
CONDA_CMD <- Sys.getenv("IMIC_CONDA_CMD", "C:/Users/andre/miniconda3/Scripts/conda.exe"); CONDA_ENV <- "mummichog"
MZ_PPM <- 10; NET_LIB <- "human_mfn"; P_CUTOFF <- 0.05   # = MUM_PPM/MUM_NET/MUM_CUTOFF in _blood_helpers.R

# --- m/z-RT references ----------------------------------------------------------
# NOTE: these feat_meta_* extractors are kept LOCAL (not the helper mzrt_* ones) on
# purpose: the directional join needs ionization mode carried as a separate column
# derived from the result feature name, so the helper's auto `mode` column would clash.
feat_meta_rlc <- function(file) {
  d <- fread(paste0(add, file), select = c("MZ","RT","Metabolite_Feature_Label"))
  lab <- toupper(d$Metabolite_Feature_Label)
  data.table(feature = lab, mz = d$MZ, rt = d$RT,
             mode = data.table::fifelse(grepl("_POS_", lab), "positive",
                    data.table::fifelse(grepl("_NEG_", lab), "negative", NA_character_)))
}
feat_meta_vam <- function() {
  d <- fread(paste0(add, "metabolite_description_vam_with_global_id.csv"))
  data.table(feature = toupper(d$feature_label), mz = d$mz, rt = d$rt_minute, mode = tolower(d$ionization_mode))
}
align_milk <- function() {
  a <- fread(paste0(add, "IMiC_alignment.csv"))
  rbind(data.table(feature = toupper(a$mtb_id_MISAME3),            mz = a$mz_MISAME3,            rt = a$rt_MISAME3),
        data.table(feature = toupper(a$mtb_id_CHILD_ELICIT_VITAL), mz = a$mz_CHILD_ELICIT_VITAL, rt = a$rt_CHILD_ELICIT_VITAL)
  )[!is.na(mz) & !is.na(rt) & feature != ""] |> unique()
}

# --- build feature tables (feature, est, pval, mz, rt, mode) per compartment ----
blood <- readRDS(paste0(root, "results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS"))
milk  <- readRDS(paste0(root, "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS"))
META  <- list(MaternalPlasma = feat_meta_rlc("ProcessedDataMISAME3_plasma.csv"),
              VamsPostnatalMaternal = feat_meta_vam(), VamsPostnatalInfant = feat_meta_vam())
ALIGN <- align_milk()

rep_feat <- function(df) df %>% group_by(feature) %>% slice_min(pval, n = 1, with_ties = FALSE) %>% ungroup()
build_blood <- function(ds, visits) {
  r <- blood %>% filter(dataset == ds, visit %in% visits, measure == "ATE") %>%
    mutate(feature = toupper(biomarker)) %>% rep_feat() %>% transmute(feature, est, pval)
  inner_join(r, META[[ds]], by = "feature")
}
build_milk <- function() {
  r <- milk %>% filter(measure == "ATE", study == "Misame", !is.na(pval), pval > 0) %>%
    mutate(feature = toupper(biomarker), mode = ifelse(grepl("_POS_", feature), "positive", "negative")) %>%
    rep_feat() %>% transmute(feature, est, pval, mode)
  inner_join(r, ALIGN, by = "feature")
}

# --- prep + run UP/DOWN for one compartment ------------------------------------
run_dir <- function(label, df) {
  for (ion in c("positive", "negative")) for (dir in c("up", "down")) {
    sub <- df %>% filter(mode == ion, !is.na(mz), !is.na(rt), !is.na(pval), pval > 0)
    if (nrow(sub) == 0) next
    sig <- if (dir == "up") sub$est > 0 else sub$est < 0
    p_use <- ifelse(sig, sub$pval, 1)                 # only this-direction features are "significant"
    mi <- data.table(mz = sub$mz, rt = sub$rt, pval = p_use, stat = -log10(pmax(p_use, 1e-300)))
    name <- paste0(label, "_", ion, "_", dir)
    f <- paste0(indir, "/", name, ".txt"); fwrite(mi, f, sep = "\t", col.names = FALSE)
    old <- setwd(outroot)
    system2(CONDA_CMD, c("run","-n",CONDA_ENV,"mummichog","-f",normalizePath(f),"-o",name,
                         "-m",ion,"-u",MZ_PPM,"-n",NET_LIB,"-c",P_CUTOFF), stdout = TRUE, stderr = TRUE)
    setwd(old)
    ok <- length(list.files(outroot, pattern = paste0("mcg_pathwayanalysis_", name), recursive = TRUE)) > 0
    cat(sprintf("  [%-34s] sig=%d %s\n", name, sum(sig), if (ok) "OK" else "NO OUTPUT"))
  }
}

cat("Directional mummichog:\n")
run_dir("Milk", build_milk())
run_dir("MaternalPlasma", build_blood("MaternalPlasma", "pn12"))
run_dir("MaternalVAMS",  build_blood("VamsPostnatalMaternal", "pn56"))
run_dir("InfantVAMS",    build_blood("VamsPostnatalInfant", c("pn12","pn34","pn56")))

# --- collect: pathway x compartment x direction --------------------------------
suppressMessages(library(readxl))
fs <- list.files(outroot, pattern = "mcg_pathwayanalysis_.*\\.xlsx$", recursive = TRUE, full.names = TRUE)
res <- bind_rows(lapply(fs, function(f) {
  d <- read_excel(f)
  # run folder is "<timestamp>.<compartment>_<mode>_<direction>"; strip the timestamp,
  # then split the remaining tag back into compartment and direction.
  source_tag <- sub("^[0-9.]+\\.", "", basename(dirname(dirname(f))))
  d$compartment <- sub("_(positive|negative)_(up|down)$", "", source_tag)
  d$direction <- sub(".*_(up|down)$", "\\1", source_tag)
  d %>% rename(p = `p-value`) %>% select(pathway, p, overlap_size, pathway_size, compartment, direction)
}))
# per pathway x compartment x direction: best p; flag enriched (p<0.05)
summ <- res %>% group_by(pathway, compartment, direction) %>% summarise(p = min(p), .groups = "drop") %>%
  mutate(enr = p < 0.05)
write.csv(summ, paste0(root, "results/directional_pathways.csv"), row.names = FALSE)
cat("\nWrote results/directional_pathways.csv\n")
