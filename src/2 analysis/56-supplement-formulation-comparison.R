# =============================================================================
# 56-supplement-formulation-comparison.R
#
# Reviewer-1 ammunition (Andrew x Trenton call, 2026-08): the two study
# supplements were run on the same VAMS mass-spec platform, so we can ask which
# metabolites DIFFER between the MISAME-III product (APSE) and the Mumta-LW
# product (Maamta). This is a NEW analysis distinct from 25 (BEP self-composition)
# and 54 (APSE detection tiers) -- here we contrast the two formulations.
#
# Method (per Trenton's steer on the call):
#   * Mean intensity per feature over DETECTED replicates only -- no imputation
#     (APSE: 18 replicates across two donor batches, pooled; Maamta: 9).
#   * log2 fold-change = log2(mean_Maamta / mean_APSE)  (Mumta relative to MISAME).
#   * "Different" = the ratio is a robust OUTLIER vs the bulk of features (Trenton:
#     "metabolites outside a certain standard deviation"): robust z on log2FC via
#     median / MAD, two-sided pseudo p = 2*pnorm(-|z|). Features present in only
#     one product are the clearest differences and are floored into the foreground.
#   * Mummichog (same recipe as the milk pipeline: 10 ppm, human_mfn, p<0.05) on
#     the outlier foreground vs the full supplement metabolome background, per ion
#     mode, to name the pathways that distinguish the two formulations.
#
# CAVEAT (reported, not hidden): these are TECHNICAL replicates of one physical
# product each, so a per-feature t-test would conflate instrument precision with a
# real formulation difference and flag almost everything. We therefore rank by
# effect size (fold-change magnitude), not a replicate t-test.
#
# Out: results/supplement_formulation_comparison.csv       (per-feature contrast)
#      results/supplement_formulation_mummichog_pathways.csv (distinguishing pathways)
#      results/mummichog_input_supplement_compare/          (mummichog inputs)
#      results/mummichog_output_supplement_compare/         (mummichog outputs)
# =============================================================================
suppressMessages({ library(data.table) })
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))
vams <- paste0(root, "data/additional datasets/ProcessedDataMISAME3_VAMS.csv")

indir  <- paste0(root, "results/mummichog_input_supplement_compare");  dir.create(indir,  showWarnings = FALSE, recursive = TRUE)
outroot<- paste0(root, "results/mummichog_output_supplement_compare"); dir.create(outroot, showWarnings = FALSE, recursive = TRUE)
CONDA_CMD <- Sys.getenv("IMIC_CONDA_CMD", "C:/Users/andre/miniconda3/Scripts/conda.exe"); CONDA_ENV <- "mummichog"

L2FC_OUTLIER_Z <- 2      # |robust z| >= 2  => "outside a standard deviation" (foreground)
ONLY_ONE_FLOOR_P <- 1e-4 # pseudo-p for present-in-one-product-only features

# --- load the 27 supplement columns (18 APSE MISAME + 9 Maamta Mumta) ----------
allcols     <- names(fread(vams, nrows = 0))
apse_cols   <- grep("^supplement;APSE",   allcols, value = TRUE)
maamta_cols <- grep("^supplement;Maamta", allcols, value = TRUE)
stopifnot(length(apse_cols) == 18, length(maamta_cols) == 9)

d <- fread(vams, select = c("MZ","RT","Metabolite_Feature_Label", apse_cols, maamta_cols))
setnames(d, c("MZ","RT","Metabolite_Feature_Label"), c("mz","rt","feature"))
d[, mode := ion_of(feature)]

# detected = non-NA & > 0; mean over detected only (no imputation)
A <- as.matrix(d[, ..apse_cols]);   A[!is.finite(A) | A <= 0] <- NA
M <- as.matrix(d[, ..maamta_cols]); M[!is.finite(M) | M <= 0] <- NA
d[, `:=`(n_apse = rowSums(!is.na(A)), n_maamta = rowSums(!is.na(M)),
         mean_apse = rowMeans(A, na.rm = TRUE), mean_maamta = rowMeans(M, na.rm = TRUE))]
d[is.nan(mean_apse), mean_apse := NA_real_]; d[is.nan(mean_maamta), mean_maamta := NA_real_]

# --- contrast: log2FC (Mumta/MISAME) + robust-outlier pseudo p -----------------
d[n_apse > 0 & n_maamta > 0, l2fc := log2(mean_maamta / mean_apse)]
center <- median(d$l2fc, na.rm = TRUE)
scale_ <- 1.4826 * median(abs(d$l2fc - center), na.rm = TRUE)   # robust SD (MAD)
d[, z := (l2fc - center) / scale_]
d[, p := 2 * pnorm(-abs(z))]

# present-in-one-product-only: clearest differences -> floor into the foreground,
# signed toward whichever product detects it.
d[n_apse > 0 & n_maamta == 0, `:=`(z = -Inf, p = ONLY_ONE_FLOOR_P)]   # MISAME only
d[n_apse == 0 & n_maamta > 0, `:=`(z =  Inf, p = ONLY_ONE_FLOOR_P)]   # Mumta only

d[, difference := fcase(
  n_apse > 0 & n_maamta == 0,                         "MISAME only",
  n_apse == 0 & n_maamta > 0,                         "Mumta only",
  !is.na(z) & z >=  L2FC_OUTLIER_Z,                   "Higher in Mumta",
  !is.na(z) & z <= -L2FC_OUTLIER_Z,                   "Higher in MISAME",
  default = "Similar")]
d[, is_different := difference != "Similar"]

fwrite(d[, .(feature, mz, rt, mode, n_apse, n_maamta,
             mean_apse, mean_maamta, l2fc, robust_z = z, pseudo_p = p, difference)],
       paste0(root, "results/supplement_formulation_comparison.csv"))
cat("wrote results/supplement_formulation_comparison.csv |", nrow(d), "features\n")
cat("difference breakdown:\n"); print(d[, .N, by = difference][order(-N)])

# --- Mummichog: which pathways distinguish the two formulations ----------------
run_mode <- function(ion) {
  sub <- d[mode == ion & is.finite(mz) & is.finite(rt) & !is.na(p)]
  if (!nrow(sub)) return(invisible())
  mi <- data.table(mz = sub$mz, rt = sub$rt, p = sub$p,
                   stat = sign(ifelse(is.finite(sub$z), sub$z, 0)) * -log10(pmax(sub$p, 1e-300)))
  f <- paste0(indir, "/supp_compare_", ion, ".txt"); fwrite(mi, f, sep = "\t", col.names = FALSE)
  old <- setwd(outroot)
  ok <- tryCatch({ system2(CONDA_CMD, c("run","-n",CONDA_ENV,"mummichog","-f", normalizePath(f),
                    "-o", paste0("supp_compare_", ion), "-m", ion, "-u", MUM_PPM, "-n", MUM_NET,
                    "-c", MUM_CUTOFF), stdout = TRUE, stderr = TRUE); TRUE },
                 error = function(e) { cat("  mummichog failed:", conditionMessage(e), "\n"); FALSE })
  setwd(old)
  cat(sprintf("  [%s] %d features, %d 'different' (p<%.2f)\n", ion, nrow(sub),
              sum(sub$p < MUM_CUTOFF, na.rm = TRUE), MUM_CUTOFF))
  invisible(ok)
}
cat("\nMummichog (formulation-distinguishing pathways):\n")
run_mode("positive"); run_mode("negative")

suppressMessages(library(readxl))
fs <- list.files(outroot, pattern = "mcg_pathwayanalysis_.*\\.xlsx$", recursive = TRUE, full.names = TRUE)
if (length(fs)) {
  res <- rbindlist(lapply(fs, function(f) {
    x <- read_excel(f); x$mode <- if (grepl("positive", f)) "positive" else "negative"
    as.data.table(x)[, .(pathway, p = `p-value`, overlap_size, pathway_size, mode)] }), fill = TRUE)
  res <- res[order(p)]
  fwrite(res, paste0(root, "results/supplement_formulation_mummichog_pathways.csv"))
  cat("\nTop formulation-distinguishing pathways:\n")
  print(head(unique(res[, .(pathway, p = signif(p,2), mode)]), 15))
  cat("\nSaved: results/supplement_formulation_mummichog_pathways.csv\n")
} else cat("\nNo mummichog output produced -- check conda env.\n")
