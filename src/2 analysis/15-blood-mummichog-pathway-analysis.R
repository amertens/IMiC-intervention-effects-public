# =============================================================================
# 15-blood-mummichog-pathway-analysis.R
#
# Mummichog pathway enrichment for the blood-compartment metabolome, adapted from
# Trenton's milk workflow (trenton scripts/imicUntargetedMetabolomicsMummichog.Rmd)
# so the blood analysis matches the manuscript's milk method exactly.
#
# Trenton's recipe (replicated here for blood):
#   * inputs SPLIT BY ION MODE, tab-delimited, NO header, 4 cols: mz  rt  pval  stat
#     where stat = -log10(pval); ALL ATE features fed (mummichog -c defines the
#     significant set vs the full-feature background).
#   * call: conda run -n <env> mummichog -f in -o out -m <mode> -u 10 -n human_mfn -c 0.05
#   * read tables/mcg_pathwayanalysis_*.xlsx, enrichment_ratio = overlap/pathway_size,
#     then FDR across pathways within each compartment x mode.
#
# Blood m/z-RT (vs Trenton's IMiC_alignment for milk): plasma/prenatal-VAMS from
# ProcessedDataMISAME3_*.csv; postnatal VAMS from metabolite_description_vam_*.csv.
# =============================================================================

suppressMessages({library(dplyr); library(data.table)})
root <- paste0(here::here(), "/"); add <- paste0(root, "data/additional datasets/")
if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- FALSE
.osuf <- if (BLOOD_ADJUST) "_adjusted" else ""    # adjusted run uses separate input/output dirs
indir  <- paste0(root, "results/mummichog_input",  .osuf); dir.create(indir,  showWarnings = FALSE, recursive = TRUE)
outroot<- paste0(root, "results/mummichog_output", .osuf); dir.create(outroot, showWarnings = FALSE, recursive = TRUE)

# ---- Mummichog run config (match Trenton; edit for this machine) ---------------
# mummichog is a python package -> two ways to call it:
#   (A) conda (Trenton): set USE_CONDA<-TRUE, CONDA_ENV to the env holding mummichog.
#   (B) direct (pip install mummichog): set USE_CONDA<-FALSE; `mummichog` must be on PATH.
# (override-aware: a wrapper can preset these before sourcing, to keep machine-
#  specific paths out of the committed defaults.)
if (!exists("RUN_MUMMICHOG"))    RUN_MUMMICHOG <- FALSE     # TRUE once mummichog is installed
if (!exists("MUMMICHOG_ENGINE")) MUMMICHOG_ENGINE <- "cli"  # "cli" = mummichog CLI; "metaboanalystr" = R fallback
if (!exists("USE_CONDA"))        USE_CONDA <- TRUE
if (!exists("CONDA_CMD"))        CONDA_CMD <- "conda"       # full path if not on PATH (e.g. ~/miniconda3/Scripts/conda.exe)
if (!exists("CONDA_ENV"))        CONDA_ENV <- "mummichog"   # conda env holding the `mummichog` python package
if (!exists("MUMMICHOG"))        MUMMICHOG <- "mummichog"   # direct command (USE_CONDA=FALSE)
MZ_PPM   <- 10                         # -u
NET_LIB  <- "human_mfn"                # -n  (human metabolic network, as in the milk paper)
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

# ---- write Trenton-format Mummichog inputs for one compartment x visit ----------
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
      transmute(mz, rt, pval, stat)                 # 4 cols, exact order Trenton uses
    if (nrow(mi) == 0) next
    f <- paste0(indir, "/", comp, "_", vis, "_", ion, ".txt")
    fwrite(mi, f, sep = "\t", col.names = FALSE)    # NO header (Trenton: CRITICAL)
    cat(sprintf("  %-26s %-5s %-9s n=%5d  sig(p<%.2f)=%d -> %s\n",
                comp, vis, ion, nrow(mi), P_CUTOFF, sum(mi$pval < P_CUTOFF), basename(f)))
    paths <- c(paths, f)
  }
  paths
}

# ---- ENGINE A: Trenton's mummichog CLI (mode inferred from filename) ------------
run_mummichog_cli <- function(input_file) {
  mode <- if (grepl("_positive\\.txt$", input_file)) "positive" else "negative"
  name <- sub("\\.txt$", "", basename(input_file))   # mummichog uses -o as a folder NAME in CWD
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

# ---- ENGINE B: MetaboAnalystR PerformPSEA (R fallback, no python needed) --------
# Same params as the CLI (10 ppm, hsa_mfn, p<0.05, per ion mode). MetaboAnalystR
# reimplements mummichog ("mum","v2"); results are very close but not identical to
# the v1 CLI used for the milk paper -- prefer the CLI for strict method-matching.
# NOTE: verify SetPeakFormat()/library codes against your MetaboAnalystR version.
run_mummichog_metaboanalyst <- function(input_file) {
  suppressMessages(library(MetaboAnalystR))
  mode <- if (grepl("_positive\\.txt$", input_file)) "positive" else "negative"
  # re-emit our 4-col file with the header MetaboAnalystR expects (m/z, rt, p, t)
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

# ---- collect pathway tables (Trenton reads tables/mcg_pathwayanalysis_*.xlsx) ---
collect_results <- function() {
  fs <- list.files(outroot, pattern = "mcg_pathwayanalysis_.*\\.xlsx$", recursive = TRUE, full.names = TRUE)
  if (length(fs) == 0) { message("No mummichog pathway tables yet."); return(invisible(NULL)) }
  # mummichog writes a NEW <timestamp>.<source> folder per run and never overwrites, so a
  # re-run leaves stale folders that would double pathway rows and corrupt the per-source
  # FDR. Keep only the LATEST-timestamp folder per source.
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
cat("Writing Mummichog inputs (Trenton format):\n")
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
