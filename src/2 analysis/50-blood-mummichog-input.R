# =============================================================================
# 50-blood-mummichog-input.R
#
# Prepares the per-feature blood effect tables used for Mummichog annotation and
# for the supplement-detection step in script 54. For three MISAME-III compartments
# (maternal plasma, maternal postnatal VAMS, infant postnatal VAMS; covariate-
# adjusted combined arms) it lists every feature at every visit, with no
# significance filter because Mummichog needs the full list as its background:
# time point, contrast, effect size (ATE), direction, raw p-value, the signed
# statistic sign(est) * -log10(p), and the m/z, retention time and ionization mode
# from each compartment's feature catalogue. Script 54 recomputes the per
# compartment x visit FDR from this table and checks the significant-feature counts.
#
# Inputs : results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS
#          data/additional datasets/{ProcessedDataMISAME3_plasma.csv,
#          metabolite_description_vam_with_global_id.csv} (via _blood_helpers.R)
# Outputs: results/blood_mummichog_input.csv    (one tidy table; read by script 54)
#          results/blood_mummichog_input.xlsx   (one sheet per compartment)
#          results/mummichog_input_blood/<comp>_<visit>_<mode>.txt
#            (tab-delimited peak lists with a header row: mz, rt, p.value, t.score)
# [needs restricted data]
# =============================================================================
suppressMessages({ library(data.table); library(openxlsx) })
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))

UP <- function(x) toupper(as.character(x))

# --- inputs -----------------------------------------------------------------
# Covariate-adjusted combined-arm blood ATEs (pooled postnatal-BEP vs control), as
# used for the cross-compartment results.
blood <- as.data.table(readRDS(paste0(
  root, "results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))
blood <- blood[measure == "ATE" & !is.na(pval) & pval > 0 & !is.na(est)]

# Each compartment's m/z-RT catalogue. Maternal + infant postnatal VAMS share the
# one V3 `vam_` catalogue; maternal plasma is its own rLC table.
CATALOGUE <- list(
  MaternalPlasma        = function() mzrt_rlc("ProcessedDataMISAME3_plasma.csv"),
  VamsPostnatalMaternal = function() mzrt_vams(),
  VamsPostnatalInfant   = function() mzrt_vams())

COMP_LABEL <- c(MaternalPlasma = "Maternal plasma",
                VamsPostnatalMaternal = "Maternal blood (postnatal VAMS)",
                VamsPostnatalInfant = "Infant blood (postnatal VAMS)")

# --- build one tidy per-feature table --------------------------------------
build_one <- function(comp) {
  d <- blood[dataset == comp]
  if (!nrow(d)) return(NULL)
  key <- CATALOGUE[[comp]]()                        # feature, mz, rt, (mode)
  setDT(key); key[, feature := UP(feature)]
  # VAMS carries ionization mode as a `mode` column (the id has no _POS_/_NEG_);
  # rLC ids encode it in the name. Keep the catalogue mode when present.
  has_mode <- "mode" %in% names(key)
  key <- unique(key[, c("feature", "mz", "rt", if (has_mode) "mode"), with = FALSE])

  d[, feature := UP(biomarker)]
  d <- merge(d, key, by = "feature", all.x = TRUE)
  d[, ion_mode := if (has_mode) fifelse(!is.na(mode), as.character(mode), ion_of(feature))
                  else ion_of(feature)]
  d[, `:=`(
    compartment = COMP_LABEL[comp],
    comp_key    = comp,
    direction   = ifelse(est < 0, "down", "up"),
    tscore      = signed_stat(pval, est))]           # sign(est) * -log10(p): the Mummichog stat

  d[, .(compartment, comp_key, timepoint = visit, contrast, feature,
        mz, rt, ion_mode, effect_size = est, direction, raw_p = pval, tscore)][order(timepoint, raw_p)]
}

tabs <- lapply(names(CATALOGUE), build_one)
names(tabs) <- names(CATALOGUE)
tidy <- rbindlist(tabs, use.names = TRUE)

# comp_key is internal (used only for the peak-file names); drop it from the tables.
deliver_cols <- setdiff(names(tidy), "comp_key")

# --- write the tidy table + a workbook with one sheet per compartment -------
outdir <- paste0(root, "results")
fwrite(tidy[, ..deliver_cols], file.path(outdir, "blood_mummichog_input.csv"))

wb <- createWorkbook()
for (comp in names(tabs)) {
  if (is.null(tabs[[comp]])) next
  sheet <- substr(comp, 1, 31)
  addWorksheet(wb, sheet); writeData(wb, sheet, tabs[[comp]][, ..deliver_cols])
}
saveWorkbook(wb, file.path(outdir, "blood_mummichog_input.xlsx"), overwrite = TRUE)

# --- write Mummichog peak lists (mz, rt, p.value, t.score) ------------------
# One tab-delimited file per compartment x visit x ionization mode, with a header
# row; only features with a mass and retention time.
mcdir <- file.path(outdir, "mummichog_input_blood")
dir.create(mcdir, showWarnings = FALSE, recursive = TRUE)
peak <- tidy[!is.na(mz) & !is.na(rt) & ion_mode %in% c("positive", "negative")]
for (ck in unique(peak$comp_key)) for (vis in unique(peak[comp_key == ck]$timepoint))
  for (mode in c("positive", "negative")) {
    sub <- peak[comp_key == ck & timepoint == vis & ion_mode == mode]
    if (!nrow(sub)) next
    fwrite(sub[, .(mz, rt, `p.value` = raw_p, `t.score` = tscore)],
           file.path(mcdir, paste0(ck, "_", vis, "_", mode, ".txt")), sep = "\t")
  }

# --- console summary --------------------------------------------------------
cat("\nblood Mummichog input: feature counts (all features per compartment x visit):\n")
print(tidy[, .(n_features = .N,
               with_mz = sum(!is.na(mz)),
               n_up = sum(direction == "up"), n_down = sum(direction == "down")),
           by = .(compartment, timepoint)][order(compartment, timepoint)])
cat("\nwrote:\n  results/blood_mummichog_input.csv (", nrow(tidy), " rows)\n",
    "  results/blood_mummichog_input.xlsx (", length(tabs), " sheets)\n",
    "  results/mummichog_input_blood/*.txt (peak lists)\n", sep = "")
