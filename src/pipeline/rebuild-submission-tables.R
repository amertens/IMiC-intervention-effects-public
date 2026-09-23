# =============================================================================
# rebuild-submission-tables.R
#
# Refresh the sheets of Manuscript/submission_tables_S1-S11.xlsx that HAVE a
# generating artifact, in place, leaving every other sheet byte-untouched.
#
# WHY ONLY SOME SHEETS
#   The workbook has never had a builder -- it is assembled by hand, which is why
#   five of its sheets (S1, S3, S4, S5, S6) lost their header rows and read back as
#   columns "0".."9". Rewriting those from their artifacts would silently change
#   column sets and row counts in sheets this audit did not verify, including the
#   two deliberate web-tool exports (S4, S6) that must NOT be overwritten from the
#   scripted pipeline. So this script touches only the sheets whose artifact was
#   verified here:
#       S10 <- results/tables/table_s10_temporal_persistence.csv
#       S11 <- results/tables/table_s11_compartment_pathway.csv
#   Extending it to another sheet means verifying that sheet's artifact first.
#
#   2026-08-26: S5 added. Verified against results/metaboanalyst/mummichog_s5/
#   milk_mummichog_tableS5.csv (P<0.05, sorted by P ascending, 505 rows) -- the
#   same filter/sort convention used to build the printed qmd Table S5 block.
#   Column set/order matches the qmd table exactly (Pathway, Overlap Size,
#   Pathway Size, Enrichment Ratio, P-value, FDR, Study, Time Point, Contrast,
#   Regulation, Ionization Mode); pathway names are left as plain text (no
#   pandoc "~9~" subscript markup) since this is a spreadsheet, not markdown.
#
#   2026-08-26: S2 and S3 added, both regenerated after the printed qmd tables
#   were found to have drifted from current data (Table S2 was a stale web
#   export disagreeing with Fig 5B; Table S3's MISAME-III rows did not
#   reproduce -- up to 272x P-value drift -- and 16 rows had alpha-/
#   gamma-linolenic acid names swapped). Both sheets are now built from the
#   same source used for the printed qmd tables:
#     S2 <- results/metaboanalyst/tertiary_msea/tertiary_msea_dual.csv
#           (all 463 rows, unfiltered, sorted by P -- matches Table S1's
#           convention of reporting the full pathway x cell result set)
#     S3 <- results/metaboanalyst/triglyceride_fa/triglyceride_fa_composition_combined.csv
#           (ELICIT + Mumta-LW rows) UNION
#           results/metaboanalyst/triglyceride_fa/triglyceride_fa_composition_stratified.csv
#           (MISAME-III rows only -- the printed table has always mixed
#           combined-arm ELICIT/Mumta-LW with stratified-arm MISAME-III; both
#           files are produced by figure-scripts/manuscript_figures/
#           fig5C-triglyceride.R, run with IMIC_TG_ARM_FRAMING=combined and
#           =stratified respectively), filtered to P<0.05, sorted by P
#           ascending (44 rows)
#
# A timestamped backup is written next to the workbook before it is modified.
#
# Run from repo root: Rscript src/pipeline/rebuild-submission-tables.R
# =============================================================================
suppressMessages({ library(openxlsx); library(data.table) })
root <- paste0(here::here(), "/")
XLSX <- paste0(root, "Manuscript/submission_tables_S1-S11.xlsx")
stopifnot(file.exists(XLSX))

# ---- S10: printed column set, em dash for undefined cells -------------------
s10 <- fread(paste0(root, "results/tables/table_s10_temporal_persistence.csv"))
dash <- function(x) ifelse(is.na(x), "—", as.character(x))
pct  <- function(x) ifelse(is.na(x), "—",
                           ifelse(x %% 1 == 0, sprintf("%.0f", x), sprintf("%.1f", x)))
S10 <- data.frame(
  `Compartment`        = s10$compartment,
  `Visit transition`   = ifelse(is.na(s10$visit_from),
                                paste0("(baseline) → ", s10$visit_to),
                                paste0(s10$visit_from, " → ", s10$visit_to)),
  `Sig. (preceding)`   = dash(s10$n_sig_preceding),
  `Sig. (current)`     = s10$n_sig_current,
  `Retained`           = s10$n_retained,
  `Newly sig.`         = s10$n_newly_sig,
  `Retained (%)`       = pct(s10$pct_retained),
  check.names = FALSE, stringsAsFactors = FALSE)

# ---- S11: printed column set ------------------------------------------------
s11 <- fread(paste0(root, "results/tables/table_s11_compartment_pathway.csv"))
S11 <- data.frame(
  `Pathway`  = s11$pathway,
  `Total`    = s11$total,
  `Expected` = round(s11$expected, 3),
  `Hits`     = s11$hits,
  `P`        = round(s11$raw_p, 3),
  `FDR`      = round(s11$fdr, 3),
  `Impact`   = round(s11$impact, 3),
  `Driving putative metabolites` =
    ifelse(is.na(s11$driving_putative_metabolites), "—", s11$driving_putative_metabolites),
  check.names = FALSE, stringsAsFactors = FALSE)

# ---- S5: printed column set (P<0.05, sorted by P ascending) -----------------
s5 <- fread(paste0(root, "results/metaboanalyst/mummichog_s5/milk_mummichog_tableS5.csv"))
s5 <- s5[s5[["P-value"]] < 0.05][order(`P-value`)]
fmt_time <- function(x) {
  x <- gsub("^([0-9.]+)-([0-9.]+) mo\\.$", "\\1\u2013\\2 months", x)
  x <- gsub("^([0-9.]+) mo\\.$", "\\1 months", x)
  x <- gsub("^1 months$", "1 month", x)
  x <- gsub("^([0-9]+)-([0-9]+) days$", "\\1\u2013\\2 days", x)
  x
}
S5 <- data.frame(
  `Pathway`             = s5$Pathway,
  `Overlap Size`         = s5$`Overlap Size`,
  `Pathway Size`         = s5$`Pathway Size`,
  `Enrichment Ratio`     = round(s5$`Enrichment Ratio`, 3),
  `P-value`              = formatC(s5$`P-value`, format = "e", digits = 3),
  `FDR`                  = formatC(s5$FDR, format = "e", digits = 3),
  `Study`                = s5$Study,
  `Time Point`           = fmt_time(s5$`Time Point`),
  `Contrast`             = s5$Contrast,
  `Regulation`           = s5$Regulation,
  `Ionization Mode`      = s5$`Ionization Mode`,
  check.names = FALSE, stringsAsFactors = FALSE)

# ---- S2: full pathway x cell result set (unfiltered, sorted by P) -----------
s2 <- fread(paste0(root, "results/metaboanalyst/tertiary_msea/tertiary_msea_dual.csv"))
s2 <- s2[order(raw_p)]
S2 <- data.frame(
  `Pathway`           = s2$pathway,
  `Total`             = s2$total,
  `Expected`          = formatC(s2$expected, format = "e", digits = 3),
  `Hits`              = s2$hits,
  `P-value`           = formatC(s2$raw_p, format = "e", digits = 3),
  `FDR`               = formatC(s2$fdr_native, format = "e", digits = 3),
  `Enrichment Ratio`  = round(s2$enrichment_ratio, 3),
  `Study`             = s2$study,
  `Contrast`          = s2$contrast,
  `Time Point`        = s2$timepoint,
  check.names = FALSE, stringsAsFactors = FALSE)

# ---- S3: combined-arm (ELICIT/Mumta-LW) UNION stratified-arm (MISAME-III) ---
s3_comb  <- fread(paste0(root, "results/metaboanalyst/triglyceride_fa/triglyceride_fa_composition_combined.csv"))
s3_strat <- fread(paste0(root, "results/metaboanalyst/triglyceride_fa/triglyceride_fa_composition_stratified.csv"))
s3 <- rbind(s3_comb[study %in% c("Elicit", "Vital")], s3_strat[study == "Misame"])
s3[, study_label := c(Elicit = "ELICIT", Vital = "Mumta-LW", Misame = "MISAME-III")[study]]
s3 <- s3[raw_p < 0.05][order(raw_p)]
s3[is.na(fatty_acid_name), fatty_acid_name := ""]
S3 <- data.frame(
  `Fatty Acid`         = s3$fatty_acid,
  `Fatty Acid Name`    = s3$fatty_acid_name,
  `Enrichment Ratio`   = round(s3$enrichment_ratio, 3),
  `P-value`            = formatC(s3$raw_p, format = "e", digits = 3),
  `FDR`                = formatC(s3$fdr_native, format = "e", digits = 3),
  `Study`              = s3$study_label,
  `Time Point`         = s3$timepoint,
  `Regulation`         = ifelse(s3$direction == "up", "Upregulated", "Downregulated"),
  `Contrast`           = s3$contrast,
  check.names = FALSE, stringsAsFactors = FALSE)

NEW <- list(S2 = S2, S3 = S3, S5 = S5, S10 = S10, S11 = S11)

backup <- sub("\\.xlsx$", paste0("_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"), XLSX)
file.copy(XLSX, backup, overwrite = FALSE)

wb  <- loadWorkbook(XLSX)
ord <- sheets(wb)                      # preserve the existing sheet order
for (nm in names(NEW)) {
  if (!nm %in% ord) { warning("sheet ", nm, " absent; skipped"); next }
  pos <- match(nm, ord)
  # Clear and overwrite IN PLACE rather than removeWorksheet()/addWorksheet(): the
  # remove/add pair appends the sheet at the end and openxlsx's worksheetOrder<-
  # would not reliably restore the original tab layout. Both rebuilt tables are at
  # least as large as the versions they replace, so no stale cells survive; the
  # generous clear window covers them regardless.
  deleteData(wb, nm, cols = 1:40, rows = 1:max(200, nrow(NEW[[nm]]) + 10), gridExpand = TRUE)
  writeData(wb, nm, NEW[[nm]], headerStyle = createStyle(textDecoration = "bold"))
  setColWidths(wb, nm, cols = seq_len(ncol(NEW[[nm]])), widths = "auto")
  worksheetOrder(wb) <- append(setdiff(seq_along(sheets(wb)), length(sheets(wb))),
                               length(sheets(wb)), after = pos - 1)
  cat("refreshed sheet", nm, "-", nrow(NEW[[nm]]), "rows x", ncol(NEW[[nm]]), "cols\n")
}
saveWorkbook(wb, XLSX, overwrite = TRUE)
cat("wrote", XLSX, "\n  backup:", backup, "\n  sheet order:", paste(sheets(wb), collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# Online-Appendix data files for the three large tables (S2, S5, S6) that are
# externalized from supplement_v2.qmd per Editor SS3.14 -- the qmd keeps only a
# short caption + pointer, matching how Table S7 (634 rows) is already handled.
# These CSVs go into the Zenodo Online Appendix snapshot.
# ---------------------------------------------------------------------------
appdir <- paste0(root, "results/tables/")
# S2: full tertiary MSEA (all 463 rows, sorted by raw P) -- same source as the sheet
s2raw <- fread(paste0(root, "results/metaboanalyst/tertiary_msea/tertiary_msea_dual.csv"))
fwrite(s2raw[order(raw_p)], paste0(appdir, "table_s3_tertiary_msea_full.csv"))
# S5: milk mummichog, raw P < 0.05, sorted by P (505 rows) -- same filter as the printed block
s5raw <- fread(paste0(root, "results/metaboanalyst/mummichog_s5/milk_mummichog_tableS5.csv"))
fwrite(s5raw[`P-value` < 0.05][order(`P-value`)], paste0(appdir, "table_s6_mummichog_full.csv"))
# S6: proteomics GO enrichment -- the current scripted source of truth (run-proteomics-go.R)
s6raw <- fread(paste0(root, "results/metaboanalyst/proteomics_go/proteomics_go_tableS6.csv"))
fwrite(s6raw, paste0(appdir, "table_s7_proteomics_go_full.csv"))
cat("wrote Online-Appendix CSVs: table_s3_tertiary_msea_full.csv (", nrow(s2raw), " rows), ",
    "table_s6_mummichog_full.csv (", nrow(s5raw[`P-value` < 0.05]), " rows), ",
    "table_s7_proteomics_go_full.csv (", nrow(s6raw), " rows)\n", sep = "")
