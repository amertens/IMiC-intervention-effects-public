# =============================================================================
# 15-blood-mummichog-pathway-analysis.R
#
# Mummichog pathway enrichment for the MISAME-III blood metabolome (combined-arm
# ATEs), run with the same settings as the milk Mummichog analysis of the original
# R Markdown workflow. Compartment x visit runs: maternal plasma 1-2 mo, maternal
# postnatal VAMS 5-6 mo, infant VAMS 1-2, 3-4 and 5-6 mo. When run with BLOOD_ADJUST
# = TRUE (as run_blood_adjusted_downstream.R does), its run folders in
# results/mummichog_output_adjusted/ are read by script 23 (per-feature putative
# annotations) and script 32 (binomial sign test for pathway direction, Methods).
#
# Method, per compartment x visit x ionization mode:
#   * input: tab-delimited, no header, 4 columns mz, rt, pval, stat with
#     stat = -log10(pval); all ATE features are supplied (mummichog -c 0.05 defines
#     the significant set against the full-feature background).
#   * call: conda run -n <env> mummichog -f in -o out -m <mode> -u 10 -n human_mfn -c 0.05
#   * collect tables/mcg_pathwayanalysis_*.xlsx, enrichment_ratio = overlap/pathway_size,
#     then BH FDR across pathways within each run.
# m/z and RT: plasma and prenatal VAMS from ProcessedDataMISAME3_*.csv; postnatal
# VAMS from metabolite_description_vam_with_global_id.csv.
#
# Inputs : results/blood_compartment_[adjusted_]combined_arms_intervention_effects_results_clean.RDS
#          data/additional datasets/{ProcessedDataMISAME3_plasma.csv, ProcessedDataMISAME3_VAMS.csv,
#          metabolite_description_vam_with_global_id.csv}
# Outputs: results/mummichog_input[_adjusted]/, results/mummichog_output[_adjusted]/,
#          results/blood_mummichog_pathways[_adjusted].{csv,RDS}
# Requires the mummichog Python package (conda environment "mummichog"; the conda
# executable is taken from IMIC_CONDA_CMD) when RUN_MUMMICHOG = TRUE.
# [needs restricted data]
# =============================================================================

suppressMessages({library(dplyr); library(data.table)})
root <- paste0(here::here(), "/"); add <- paste0(root, "data/additional datasets/")
if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- FALSE
.osuf <- if (BLOOD_ADJUST) "_adjusted" else ""    # adjusted run uses separate input/output dirs
indir  <- paste0(root, "results/mummichog_input",  .osuf); dir.create(indir,  showWarnings = FALSE, recursive = TRUE)
outroot<- paste0(root, "results/mummichog_output", .osuf); dir.create(outroot, showWarnings = FALSE, recursive = TRUE)

# ---- Mummichog run settings ----------------------------------------------------
# mummichog is a Python package; two ways to call it:
#   (A) conda: USE_CONDA <- TRUE and CONDA_ENV = the environment holding mummichog.
#   (B) direct (pip install mummichog): USE_CONDA <- FALSE; `mummichog` must be on PATH.
# A wrapper can preset any of these before sourcing.
if (!exists("RUN_MUMMICHOG"))    RUN_MUMMICHOG <- FALSE     # TRUE once mummichog is installed
if (!exists("MUMMICHOG_ENGINE")) MUMMICHOG_ENGINE <- "cli"  # "cli" = mummichog CLI; "metaboanalystr" = R fallback
if (!exists("USE_CONDA"))        USE_CONDA <- TRUE
if (!exists("CONDA_CMD"))        CONDA_CMD <- Sys.getenv("IMIC_CONDA_CMD", "conda")  # full path if conda is not on PATH
if (!exists("CONDA_ENV"))        CONDA_ENV <- "mummichog"   # conda env holding the `mummichog` python package
if (!exists("MUMMICHOG"))        MUMMICHOG <- "mummichog"   # direct command (USE_CONDA=FALSE)
MZ_PPM   <- 10                         # -u
NET_LIB  <- "human_mfn"                # -n  (human metabolic network, as for milk)
P_CUTOFF <- 0.05                       # -c  (defines the significant feature set)

# ---- feature -> (mz, rt, ion mode) per compartment -----------------------------
feat_meta_rlc <- function(file) {
  d <- fread(paste0(add, file), select = c("MZ", "RT", "Metabolite_Feature_Label"))
  lab <- toupper(d$Metabolite_Feature_Label)
  data.table(feature = lab, mz = d$MZ, rt = d$RT,
             mode = data.table::fifelse(grepl("_POS_", lab), "positive",
                    data.table::fifelse(grepl("_NEG_", lab), "negative", NA_character_)))
}
feat_meta_vam <- function() {
  d <- fread(paste0(add, "metabolite_description_vam_with_global_id.csv"))
  data.table(feature = toupper(d$feature_label), mz = d$mz, rt = d$rt_minute, mode = tolower(d$ionization_mode))
}
META <- list(
  MaternalPlasma        = feat_meta_rlc("ProcessedDataMISAME3_plasma.csv"),
  VamsPrenatal          = feat_meta_rlc("ProcessedDataMISAME3_VAMS.csv"),
  VamsPostnatalMaternal = feat_meta_vam(),
  VamsPostnatalInfant   = feat_meta_vam())

.btag <- if (BLOOD_ADJUST) "adjusted_" else ""   # read adjusted_ blood results when set (BLOOD_ADJUST set above)
blood <- readRDS(paste0(root, "results/blood_compartment_", .btag,
  "combined_arms_intervention_effects_results_clean.RDS"))

# ---- write the 4-column Mummichog inputs for one compartment x visit -----------
# (params comp/vis/ion to avoid colliding with the data columns of the same name.)
prepare_input <- function(comp, vis) {
  res <- blood %>%
    filter(dataset == comp, visit == vis, measure == "ATE", !is.na(pval), pval > 0) %>%
    mutate(feature = toupper(biomarker), stat = -log10(pval)) %>%
    transmute(feature, pval, stat)
  m <- inner_join(res, META[[comp]], by = "feature")
  paths <- character()
  for (ion in c("positive", "negative")) {
    mi <- m %>% filter(mode == ion, !is.na(mz), !is.na(rt)) %>%
      transmute(mz, rt, pval, stat)                 # column order mummichog expects
    if (nrow(mi) == 0) next
    f <- paste0(indir, "/", comp, "_", vis, "_", ion, ".txt")
    fwrite(mi, f, sep = "\t", col.names = FALSE)    # no header row, as the milk inputs
    cat(sprintf("  %-26s %-5s %-9s n=%5d  sig(p<%.2f)=%d -> %s\n",
                comp, vis, ion, nrow(mi), P_CUTOFF, sum(mi$pval < P_CUTOFF), basename(f)))
    paths <- c(paths, f)
  }
  paths
}

# ---- engine A: mummichog CLI (ion mode inferred from the file name) -------------
run_mummichog_cli <- function(input_file) {
  mode <- if (grepl("_positive\\.txt$", input_file)) "positive" else "negative"
  name <- sub("\\.txt$", "", basename(input_file))   # mummichog uses -o as a folder name in the working directory
  inf  <- normalizePath(input_file)
  old  <- setwd(outroot); on.exit(setwd(old))         # output created (timestamped) relative to CWD
  mflags <- c("-f", inf, "-o", name, "-m", mode, "-u", MZ_PPM, "-n", NET_LIB, "-c", P_CUTOFF)
  cmd  <- if (USE_CONDA) CONDA_CMD else MUMMICHOG
  args <- if (USE_CONDA) c("run", "-n", CONDA_ENV, "mummichog", mflags) else mflags
  t <- Sys.time()
  invisible(tryCatch(system2(cmd, args = args, stdout = TRUE, stderr = TRUE),
                     error = function(e) paste("ERROR:", conditionMessage(e))))
  ok <- length(list.files(outroot, pattern = paste0("mcg_pathwayanalysis_", name), recursive = TRUE)) > 0
  cat(sprintf("  [%-32s] %s  (%.0fs)\n", name, if (ok) "OK" else "NO OUTPUT",
              as.numeric(difftime(Sys.time(), t, units = "secs"))))
  name
}

# ---- engine B: MetaboAnalystR PerformPSEA (R fallback, no Python needed) --------
# Same parameters as the CLI (10 ppm, hsa_mfn, p < 0.05, per ion mode). MetaboAnalystR
# reimplements mummichog ("mum", "v2"); results are close to but not identical to the
# v1 CLI used for milk, so the CLI is the engine used for the reported results.
# SetPeakFormat() and library codes can differ between MetaboAnalystR versions.
run_mummichog_metaboanalyst <- function(input_file) {
  suppressMessages(library(MetaboAnalystR))
  mode <- if (grepl("_positive\\.txt$", input_file)) "positive" else "negative"
  # re-write the 4-column file with the header MetaboAnalystR expects (m/z, rt, p, t)
  d <- data.table::fread(input_file, header = FALSE,
                         col.names = c("m.z", "r.t", "p.value", "t.score"))
  ma_in <- sub("\\.txt$", "_ma.txt", input_file)
  data.table::fwrite(d, ma_in, sep = "\t")
  out <- paste0(outroot, "/", sub("\\.txt$", "", basename(input_file)))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  old <- setwd(out); on.exit(setwd(old))
  mSet <- InitDataObjects("mass_all", "mummichog", FALSE)
  mSet <- SetPeakFormat(mSet, "rmp")                              # rt + m/z + p (+ t)
  mSet <- UpdateInstrumentParameters(mSet, MZ_PPM, mode, "yes", 0.02)
  mSet <- Read.PeakListData(mSet, ma_in)
  mSet <- SanityCheckMummichogData(mSet)
  mSet <- SetPeakEnrichMethod(mSet, "mum", "v2")
  mSet <- SetMummichogPval(mSet, P_CUTOFF)
  mSet <- PerformPSEA(mSet, "hsa_mfn", "current", 3, 100)         # hsa_mfn == CLI human_mfn
  mSet
}

# dispatcher
run_mummichog <- function(input_file) {
  if (MUMMICHOG_ENGINE == "metaboanalystr") run_mummichog_metaboanalyst(input_file)
  else run_mummichog_cli(input_file)
}

# ---- collect pathway tables (tables/mcg_pathwayanalysis_*.xlsx) -----------------
collect_results <- function() {
  fs <- list.files(outroot, pattern = "mcg_pathwayanalysis_.*\\.xlsx$", recursive = TRUE, full.names = TRUE)
  if (length(fs) == 0) { message("No mummichog pathway tables yet."); return(invisible(NULL)) }
  # mummichog writes a new <timestamp>.<source> folder per run and never overwrites, so a
  # re-run leaves stale folders that would double pathway rows and distort the per-source
  # FDR. Keep only the latest-timestamp folder per source.
  fold <- basename(dirname(dirname(fs)))                                  # "<ts>.<source>"
  pick <- data.table::data.table(
            f  = fs,
            src = sub("^[0-9]+\\.[0-9]+\\.", "", fold),
            ts  = suppressWarnings(as.numeric(sub("^([0-9]+\\.[0-9]+)\\..*$", "\\1", fold)))
          )[order(-ts)][!duplicated(src)]
  fs <- pick$f
  suppressMessages(library(readxl))
  res <- bind_rows(lapply(fs, function(f) {
    d <- readxl::read_excel(f)
    d$source <- sub("^[0-9.]+\\.", "", basename(dirname(dirname(f))))  # strip mummichog timestamp prefix
    d
  })) %>%
    mutate(enrichment_ratio = overlap_size / pathway_size) %>%
    group_by(source) %>% mutate(pathway_FDR = p.adjust(`p-value`, method = "BH")) %>% ungroup()
  saveRDS(res, paste0(root, "results/blood_mummichog_pathways", .osuf, ".RDS"))
  write.csv(res, paste0(root, "results/blood_mummichog_pathways", .osuf, ".csv"), row.names = FALSE)
  cat("Collected", nrow(res), "pathway rows from", length(fs), "runs.\n")
  invisible(res)
}

# ---- driver --------------------------------------------------------------------
cat("Writing Mummichog inputs (mz, rt, pval, stat):\n")
inputs <- c(
  prepare_input("MaternalPlasma", "pn12"),
  prepare_input("VamsPostnatalMaternal", "pn56"),
  unlist(lapply(c("pn12", "pn34", "pn56"), function(v) prepare_input("VamsPostnatalInfant", v))))

if (RUN_MUMMICHOG) {
  cat("\nRunning mummichog (", CONDA_CMD, "run -n", CONDA_ENV, "):\n")
  for (f in inputs) run_mummichog(f)
  collect_results()
} else {
  cat("\nInputs ready in results/mummichog_input/. To run: install the `mummichog` python\n",
      "package in a conda env, set RUN_MUMMICHOG<-TRUE + CONDA_CMD/CONDA_ENV, re-source.\n",
      "Then collect_results() tidies tables/mcg_pathwayanalysis_*.xlsx with FDR across pathways.\n")
}
