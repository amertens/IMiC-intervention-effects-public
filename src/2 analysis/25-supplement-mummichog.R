# =============================================================================
# 25-supplement-mummichog.R
#
# Andrew x Trenton ask ("do the same Mummichog on the BEP"): characterise the BEP
# product itself, which metabolic pathways are the supplement's own metabolites
# enriched for? The supplement is on the V1 rLC catalogue (recovered, 27 samples);
# we treat the features that are ABUNDANT in the product (high supp_mean) as the
# "significant" set and run the same Mummichog pathway enrichment used elsewhere.
# This describes the product's composition; it is NOT an intervention contrast.
#
# Out: results/supplement_mummichog_pathways.csv  (+ mummichog_output_supplement/)
# =============================================================================

suppressMessages({library(data.table)})
root <- paste0(here::here(), "/"); add <- paste0(root, "data/additional datasets/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))
indir  <- paste0(root, "results/mummichog_input_supplement");  dir.create(indir,  showWarnings = FALSE, recursive = TRUE)
outroot<- paste0(root, "results/mummichog_output_supplement"); dir.create(outroot, showWarnings = FALSE, recursive = TRUE)
CONDA_CMD <- Sys.getenv("IMIC_CONDA_CMD", "C:/Users/andre/miniconda3/Scripts/conda.exe"); CONDA_ENV <- "mummichog"

# supplement abundance (V1 rLC) + m/z-RT-mode from the prenatal-VAMS V1 catalogue
supp <- fread(paste0(root, "results/bep_supplement_profile_V1.csv"))
pv   <- mzrt_rlc("ProcessedDataMISAME3_VAMS.csv")
d <- merge(data.table(feature = toupper(supp$feature), supp_mean = supp$supp_mean),
           pv[, .(feature, mz, rt, mode)], by = "feature")
# pseudo p-value: higher product abundance -> more "significant" (one-sided normal tail)
d[, p := pnorm(-supp_mean)]
d <- d[is.finite(mz) & is.finite(rt) & !is.na(mode) & p > 0]

run_mode <- function(ion) {
  sub <- d[mode == ion]; if (!nrow(sub)) return(invisible())
  mi <- data.table(mz = sub$mz, rt = sub$rt, p = sub$p, stat = -log10(pmax(sub$p, 1e-300)))
  f <- paste0(indir, "/supplement_", ion, ".txt"); fwrite(mi, f, sep = "\t", col.names = FALSE)
  old <- setwd(outroot)
  system2(CONDA_CMD, c("run","-n",CONDA_ENV,"mummichog","-f", normalizePath(f), "-o", paste0("supplement_", ion),
                       "-m", ion, "-u", MUM_PPM, "-n", MUM_NET, "-c", MUM_CUTOFF), stdout = TRUE, stderr = TRUE)
  setwd(old)
  cat(sprintf("  [%s] %d features, %d 'significant' (p<%.2f)\n", ion, nrow(sub), sum(sub$p < MUM_CUTOFF), MUM_CUTOFF))
}
cat("Supplement Mummichog (product composition, V1 rLC):\n")
run_mode("positive"); run_mode("negative")

suppressMessages(library(readxl))
fs <- list.files(outroot, pattern = "mcg_pathwayanalysis_.*\\.xlsx$", recursive = TRUE, full.names = TRUE)
if (length(fs)) {
  res <- rbindlist(lapply(fs, function(f) {
    x <- read_excel(f); x$mode <- if (grepl("positive", f)) "positive" else "negative"
    as.data.table(x)[, .(pathway, p = `p-value`, overlap_size, pathway_size, mode)] }), fill = TRUE)
  res <- res[order(p)]
  fwrite(res, paste0(root, "results/supplement_mummichog_pathways.csv"))
  cat("\nTop supplement-composition pathways:\n"); print(head(unique(res[, .(pathway, p = signif(p,2), mode)]), 15))
  cat("\nSaved: results/supplement_mummichog_pathways.csv\n")
} else cat("\nNo mummichog output produced: check conda env.\n")
