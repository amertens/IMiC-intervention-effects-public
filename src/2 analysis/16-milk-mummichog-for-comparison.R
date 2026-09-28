# =============================================================================
# 16-milk-mummichog-for-comparison.R
#
# Runs Mummichog on the MISAME-III milk untargeted metabolome (covariate-adjusted
# combined-arm ATEs), per visit and ionization mode, with the same settings as the
# blood runs in script 15, so milk and blood pathways are comparable. Script 33 reads
# these runs for the milk binomial sign test (Methods), and script 41 tabulates the
# result by visit: the Discussion statement that de novo fatty-acid biosynthesis was
# decreased in milk. Milk m/z and RT come from the cross-study alignment key
# (mtb_id_MISAME3 -> mz_MISAME3 / rt_MISAME3).
#
# Inputs : results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS
#          data/additional datasets/IMiC_alignment.csv
# Outputs: results/mummichog_input/Milk_<visit>_<mode>.txt
#          results/mummichog_output/<timestamp>.Milk_<visit>_<mode>/ (mummichog run folders)
# Requires the mummichog Python package in the conda environment "mummichog"; the
# conda executable is taken from IMIC_CONDA_CMD (default "conda").
# [needs restricted data]
# =============================================================================

suppressMessages({library(dplyr); library(data.table)})
root <- paste0(here::here(), "/"); add <- paste0(root, "data/additional datasets/")
indir   <- paste0(root, "results/mummichog_input");  dir.create(indir, showWarnings = FALSE, recursive = TRUE)
outroot <- paste0(root, "results/mummichog_output"); dir.create(outroot, showWarnings = FALSE, recursive = TRUE)

# Mummichog settings, as in script 15: 10 ppm, human_mfn network, p < 0.05 cutoff
CONDA_CMD <- Sys.getenv("IMIC_CONDA_CMD", "conda")
CONDA_ENV <- "mummichog"; MZ_PPM <- 10; NET_LIB <- "human_mfn"; P_CUTOFF <- 0.05

# milk untargeted ATE (MISAME-III) + alignment-key m/z-RT
align <- fread(paste0(add, "IMiC_alignment.csv")) %>%
  distinct(mtb_id_MISAME3, mz_MISAME3, rt_MISAME3) %>%
  transmute(biomarker = toupper(mtb_id_MISAME3), mz = mz_MISAME3, rt = rt_MISAME3)
milk <- readRDS(paste0(root, "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")) %>%
  filter(measure == "ATE", study == "Misame", !is.na(pval), pval > 0) %>%
  mutate(biomarker = toupper(biomarker),
         ion = data.table::fifelse(grepl("_POS_", biomarker), "positive",
               data.table::fifelse(grepl("_NEG_", biomarker), "negative", NA_character_))) %>%
  inner_join(align, by = "biomarker")

# One mummichog CLI run; mummichog creates its (timestamp-prefixed) output folder in
# the working directory, so switch to outroot for the call.
run_cli <- function(input_file) {
  mode <- if (grepl("_positive\\.txt$", input_file)) "positive" else "negative"
  name <- sub("\\.txt$", "", basename(input_file)); inf <- normalizePath(input_file)
  old <- setwd(outroot); on.exit(setwd(old))
  t <- Sys.time()
  system2(CONDA_CMD, c("run", "-n", CONDA_ENV, "mummichog", "-f", inf, "-o", name,
                       "-m", mode, "-u", MZ_PPM, "-n", NET_LIB, "-c", P_CUTOFF),
          stdout = TRUE, stderr = TRUE)
  ok <- length(list.files(outroot, pattern = paste0("mcg_pathwayanalysis_", name), recursive = TRUE)) > 0
  cat(sprintf("  [%-26s] %s (%.0fs)\n", name, if (ok) "OK" else "NO OUTPUT",
              as.numeric(difftime(Sys.time(), t, units = "secs"))))
}

cat("Milk mummichog (per visit x ion mode):\n")
for (v in unique(milk$visit)) {
  vlab <- gsub("[^A-Za-z0-9]", "", v)                      # filename-safe visit ("1-2 mo." -> "12mo")
  for (mode in c("positive", "negative")) {
    mi <- milk %>% filter(visit == v, ion == mode) %>%
      transmute(mz, rt, pval, stat = -log10(pval)) %>% filter(!is.na(mz), !is.na(rt))
    if (nrow(mi) == 0) next
    f <- paste0(indir, "/Milk_", vlab, "_", mode, ".txt")
    fwrite(mi, f, sep = "\t", col.names = FALSE)
    cat(sprintf("  prep Milk_%s_%s n=%d sig=%d\n", vlab, mode, nrow(mi), sum(mi$pval < 0.05)))
    run_cli(f)
  }
}
cat("Done milk mummichog; run folders are in results/mummichog_output/ (read by script 33).\n")
