# =============================================================================
# run-milk-mummichog.R
#
# Directional Mummichog pathway analysis of the untargeted milk metabolome for
# Fig 6B and Table S6 (Table S6 lists the rows with P < 0.05). Runs are stratified
# by study x time point x contrast, by ionization mode and by direction of effect;
# Mumta-LW 1.5 months uses the combined-arm results (time point before antibiotic
# administration), all other cells the stratified-arm results. Each run calls the
# mummichog command-line tool (human_mfn network, 10 ppm, P cutoff 0.05) through
# conda; FDR is Benjamini-Hochberg across pathways within each run.
#
# Direction is handled by submitting only the in-direction features for each
# up/down run, so Mummichog's background is the in-direction feature set, as in
# the original directional runs. Setting the opposite direction to P = 1 instead
# would keep those features in the background, roughly doubling it and shifting
# every permutation P-value.
#
# Inputs : [needs restricted data]
#          results/adjusted_intervention_effects_res_untargeted_metabolomics_clean_ATE.RDS
#          results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS
#          data/additional datasets/IMiC_alignment.csv  (feature m/z and retention time)
#          conda env "mummichog" (conda command from IMIC_CONDA_CMD, default "conda")
# Outputs: results/metaboanalyst/mummichog/milk_mummichog_pathways.csv  (Fig 6B, Table S6)
#          results/mummichog_input_milk/   (per-run mummichog input files)
#          results/mummichog_output_milk/  (per-run mummichog output; reused on re-run)
#
# Run from repo root: Rscript src/metaboanalyst/run-milk-mummichog.R
# =============================================================================
suppressMessages({ library(dplyr); library(data.table) })

ROOT  <- paste0(here::here(), "/")
ADD   <- paste0(ROOT, "data/additional datasets/")
INDIR <- paste0(ROOT, "results/mummichog_input_milk")
OUTD  <- paste0(ROOT, "results/mummichog_output_milk")
dir.create(INDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTD,  showWarnings = FALSE, recursive = TRUE)

CONDA_CMD <- Sys.getenv("IMIC_CONDA_CMD", "conda")
CONDA_ENV <- "mummichog"    # conda env holding the mummichog CLI
MZ_PPM    <- 10             # m/z match tolerance in ppm (mummichog -u)
NET_LIB   <- "human_mfn"    # metabolic network model (mummichog -n)
P_CUTOFF  <- 0.05           # significance cutoff for the "hit" list (mummichog -c)

# Build the biomarker -> (m/z, retention-time) lookup for the CHILD/ELICIT/VITAL
# studies from the shared alignment table. Used by both feature loaders below.
cev_key <- function(alignment) {
  alignment %>%
    distinct(mtb_id_CHILD_ELICIT_VITAL, mz_CHILD_ELICIT_VITAL, rt_CHILD_ELICIT_VITAL) %>%
    transmute(biomarker = toupper(mtb_id_CHILD_ELICIT_VITAL),
              mz = mz_CHILD_ELICIT_VITAL, rt = rt_CHILD_ELICIT_VITAL) %>%
    filter(!is.na(biomarker), biomarker != "")
}

# Read the ionization mode out of the biomarker id (its name encodes _POS_ / _NEG_)
# and drop any feature that is tagged with neither mode.
classify_ion <- function(features) {
  features %>%
    mutate(biomarker = toupper(biomarker),
           ion = fifelse(grepl("_POS_", biomarker), "positive",
                 fifelse(grepl("_NEG_", biomarker), "negative", NA_character_))) %>%
    filter(!is.na(ion))
}

# --- untargeted ATE (stratified arms) + per-study m/z-RT alignment key ----------
load_untargeted <- function() {
  alignment <- fread(paste0(ADD, "IMiC_alignment.csv"))
  # MISAME uses its own m/z-RT alignment columns; the other studies share the CEV key.
  misame_key <- alignment %>%
    distinct(mtb_id_MISAME3, mz_MISAME3, rt_MISAME3) %>%
    transmute(biomarker = toupper(mtb_id_MISAME3), mz = mz_MISAME3, rt = rt_MISAME3) %>%
    filter(!is.na(biomarker), biomarker != "")
  cev <- cev_key(alignment)

  effects <- readRDS(paste0(ROOT, "results/adjusted_intervention_effects_res_untargeted_metabolomics_clean_ATE.RDS")) %>%
    filter(measure == "ATE", !is.na(pval), pval > 0, !is.na(est)) %>%
    classify_ion()

  # Attach m/z + RT to each feature using the study-appropriate alignment key.
  bind_rows(
    effects %>% filter(study == "Misame") %>% inner_join(misame_key, by = "biomarker"),
    effects %>% filter(study != "Misame") %>% inner_join(cev, by = "biomarker"))
}


# Some cells use the combined-arm framework rather than stratified arms: per the
# Methods, "combined-arm analyses [were] used for time points preceding antibiotic
# administration". Empirically that is Mumta-LW (Vital) 1.5 mo, contrast "BEP".
COMBINED_CELLS <- list(list(study = "Vital", visit = "1.5 mo.", contrast = "BEP"))

load_untargeted_combined <- function() {
  alignment <- fread(paste0(ADD, "IMiC_alignment.csv"))
  cev <- cev_key(alignment)
  readRDS(paste0(ROOT, "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")) %>%
    filter(measure == "ATE", !is.na(pval), pval > 0, !is.na(est), study != "Misame") %>%
    classify_ion() %>% inner_join(cev, by = "biomarker")
}

# --- run one (study, visit, contrast, ion, direction) cell ----------------------
run_cell <- function(df, study, visit, contrast, ion, dir) {
  # Features for this single cell (study x visit x contrast x ionization mode).
  cell_features <- df %>% filter(study == !!study, visit == !!visit, contrast == !!contrast,
                                 ion == !!ion, !is.na(mz), !is.na(rt))

  # Directional run: keep only the in-direction features, with their real
  # p-values, so the background is the in-direction feature set (see header).
  cell_features <- cell_features %>% filter(if (dir == "down") est < 0 else est > 0)
  if (!nrow(cell_features)) return(NULL)

  mummichog_input <- data.table(
    mz   = cell_features$mz,
    rt   = cell_features$rt,
    pval = cell_features$pval,
    stat = sign(cell_features$est) * -log10(pmax(cell_features$pval, 1e-300)))  # signed -log10 p

  run_name   <- paste0(gsub("[^A-Za-z0-9]", "", paste0(study, visit, contrast)), "_", ion, "_", dir)
  input_path <- file.path(INDIR, paste0(run_name, ".txt"))
  fwrite(mummichog_input, input_path, sep = "\t", col.names = FALSE)

  # Tag a mummichog results table with this cell's identity (kept identical
  # across the resume and fresh-run paths so downstream binding lines up).
  annotate_run <- function(dt) {
    dt %>% mutate(study = study, visit = visit, contrast = contrast,
                  ion = ion, direction = dir, n_input = nrow(mummichog_input),
                  n_sig = sum(cell_features$pval < 0.05))
  }

  # Resume: mummichog runs are expensive, so reuse a completed run if its
  # results table is already on disk (lets the grid be restarted safely).
  finished <- list.files(OUTD, pattern = paste0("^mcg_pathwayanalysis_", run_name, "\\.tsv$"),
                         recursive = TRUE, full.names = TRUE)
  if (length(finished)) return(annotate_run(fread(finished[1])))

  old <- setwd(OUTD); on.exit(setwd(old), add = TRUE)
  system2(CONDA_CMD, c("run", "-n", CONDA_ENV, "mummichog", "-f", normalizePath(input_path),
                       "-o", run_name, "-m", ion, "-u", MZ_PPM, "-n", NET_LIB, "-c", P_CUTOFF),
          stdout = TRUE, stderr = TRUE)
  results_files <- list.files(OUTD, pattern = paste0("mcg_pathwayanalysis_", run_name, ".*\\.tsv$"),
                              recursive = TRUE, full.names = TRUE)
  if (!length(results_files)) return(NULL)
  annotate_run(fread(results_files[1]))
}

# --- full grid -> pathway table -----------------------------------------------
# Enrichment ratio is overlap_size / pathway_size, signed by direction (checked
# against the published table: 31/68 = 0.456, 30/77 = 0.390). FDR is Benjamini-
# Hochberg across pathways within each run (study x visit x contrast x ion x dir).
.study_label <- function(s) c(Misame = "MISAME", Vital = "Mumta-LW", Elicit = "ELICIT")[s]

run_all_cells <- function(write = TRUE) {
  # 1. Stratified-arm cells: every (study, visit, contrast) x ion x direction.
  effects   <- load_untargeted()
  cell_grid <- effects %>% distinct(study, visit, contrast)
  pathway_tables <- list()
  for (i in seq_len(nrow(cell_grid))) for (ion in c("positive", "negative")) for (dir in c("down", "up")) {
    result <- tryCatch(
      run_cell(effects, cell_grid$study[i], cell_grid$visit[i], cell_grid$contrast[i], ion, dir),
      error = function(e) NULL)
    if (is.null(result) || !nrow(result)) next
    pathway_tables[[length(pathway_tables) + 1]] <- result
    message(sprintf("  %s %s %s %s %s -> %d pathways",
                    cell_grid$study[i], cell_grid$visit[i], cell_grid$contrast[i], ion, dir, nrow(result)))
  }

  # 2. The handful of combined-arm cells (pre-antibiotic time points).
  combined_effects <- load_untargeted_combined()
  for (cc in COMBINED_CELLS) for (ion in c("positive", "negative")) for (dir in c("down", "up")) {
    result <- tryCatch(run_cell(combined_effects, cc$study, cc$visit, cc$contrast, ion, dir),
                       error = function(e) NULL)
    if (is.null(result) || !nrow(result)) next
    pathway_tables[[length(pathway_tables) + 1]] <- result
    message(sprintf("  [combined-arm] %s %s %s %s %s -> %d pathways",
                    cc$study, cc$visit, cc$contrast, ion, dir, nrow(result)))
  }
  if (!length(pathway_tables)) return(invisible(NULL))

  # 3. Pool all cells, FDR-correct within each run, then format the table.
  pathway_table <- bind_rows(pathway_tables) %>%
    rename(p_value = `p-value`) %>%
    group_by(study, visit, contrast, ion, direction) %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>% ungroup() %>%
    mutate(enrichment_ratio = ifelse(direction == "down", -1, 1) * (overlap_size / pathway_size),
           Study = .study_label(study),
           Regulation = ifelse(direction == "down", "Downregulated", "Upregulated"),
           `Ionization Mode` = tools::toTitleCase(ion)) %>%
    transmute(Pathway = pathway, `Overlap Size` = overlap_size, `Pathway Size` = pathway_size,
              `Enrichment Ratio` = round(enrichment_ratio, 3), `P-value` = p_value, FDR = fdr,
              Study, `Time Point` = visit, Contrast = contrast, Regulation, `Ionization Mode`) %>%
    arrange(`P-value`)
  if (isTRUE(write)) {
    dir.create(paste0(ROOT, "results/metaboanalyst/mummichog"), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(pathway_table, paste0(ROOT, "results/metaboanalyst/mummichog/milk_mummichog_pathways.csv"))
  }
  pathway_table
}

if (sys.nframe() == 0) {
  pathway_table <- run_all_cells()
  message("pathway rows: ", nrow(pathway_table), " | FDR<0.05: ", sum(pathway_table$FDR < 0.05, na.rm = TRUE))
}
