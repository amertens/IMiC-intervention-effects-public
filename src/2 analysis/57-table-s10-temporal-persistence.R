# =============================================================================
# 57-table-s10-temporal-persistence.R
#
# Builds Table S10: persistence and emergence of FDR-significant up-regulated
# BEP-associated metabolite features across consecutive visits within each MISAME-III
# blood compartment (the Results statement that 100% of infant features significant
# at 1-2 months remained so at 3-4 months, and 87.5% from 3-4 to 5-6 months).
# Definitions, as in the table caption:
#   - Universe: blood bioTMLE ATEs, pooled postnatal-BEP arms vs control (combined arm
#     coding), covariate-adjusted, contrast "BEP".
#   - Significant = BH q < 0.05 within dataset x visit. The input holds only the
#     FDR-significant rows (sigFDR == 1), so a dataset x visit cell absent from it has zero.
#   - Up-regulated = est > 0.
#   - Retained = the same feature id is FDR-significantly up-regulated at both the
#     preceding and the current visit; newly significant = current - retained;
#     retention % = retained / preceding (undefined when preceding == 0).
# These are group-level feature counts over calendar time, not within-individual
# accumulation. Counts are up-regulated only: infant blood has 3 / 8 / 17 features at
# 1-2 / 3-4 / 5-6 months, against 3 / 8 / 18 in both directions (the counts checked in
# script 54), because one 5-6 month feature is down-regulated.
#
# Input  : results/blood_compartment_all_FDRsig_ATE.csv (combine_blood_results.R)
# Outputs: results/tables/table_s10_temporal_persistence.csv
#          results/tables/table_s10_temporal_persistence.md (markdown version of the table)
# [needs restricted data] (upstream blood results)
# =============================================================================
suppressMessages({ library(data.table) })
root <- paste0(here::here(), "/")

IN   <- paste0(root, "results/blood_compartment_all_FDRsig_ATE.csv")
OUTD <- paste0(root, "results/tables")
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

# Compartment -> ordered visit sequence, with the labels used in the printed table.
# Only compartments with >= 2 analysed visits in the combined-arm BEP contrast appear.
SEQ <- list(
  list(dataset = "VamsPostnatalInfant",   label = "Infant blood (VAMS)",
       visits = c(acco = "Birth", pn12 = "1-2 mo.", pn34 = "3-4 mo.", pn56 = "5-6 mo.")),
  list(dataset = "VamsPostnatalMaternal", label = "Maternal VAMS",
       visits = c(tri3 = "Third trimester", pn56 = "5-6 mo.")),
  list(dataset = "MaternalPlasma",        label = "Maternal plasma",
       visits = c(incl = "Inclusion", tri3 = "Third trimester", pn12 = "1-2 mo."))
)

d <- fread(IN)
up <- d[contrast == "BEP" & arm_coding == "combined" & adjustment == "adjusted" &
          sigFDR == 1 & est > 0]
up[, fid := tolower(biomarker)]

# feature-id set of FDR-significant up-regulated features in one (dataset, visit)
fset <- function(ds, v) unique(up[dataset == ds & visit == v, fid])

rows <- rbindlist(lapply(SEQ, function(s) {
  vs <- names(s$visits)
  rbindlist(lapply(seq_along(vs), function(i) {
    cur  <- fset(s$dataset, vs[i])
    prev <- if (i == 1) NULL else fset(s$dataset, vs[i - 1])
    ret  <- if (i == 1) 0L else length(intersect(prev, cur))
    data.table(
      compartment      = s$label,
      dataset          = s$dataset,
      visit_from       = if (i == 1) NA_character_ else unname(s$visits[i - 1]),
      visit_to         = unname(s$visits[i]),
      visit_from_code  = if (i == 1) NA_character_ else vs[i - 1],
      visit_to_code    = vs[i],
      n_sig_preceding  = if (i == 1) NA_integer_ else length(prev),
      n_sig_current    = length(cur),
      n_retained       = ret,
      n_newly_sig      = length(cur) - ret,
      pct_retained     = if (i == 1 || length(prev) == 0) NA_real_
                         else round(100 * ret / length(prev), 1))
  }))
}))

# The printed table shows a compartment's first visit as its own row only when that
# visit already has significant features (maternal plasma, Inclusion: 16). A first
# visit with none (infant Birth, maternal-VAMS third trimester) contributes only as
# the "preceding" column of the next transition, so drop those empty leading rows.
rows <- rows[!(is.na(visit_from) & n_sig_current == 0)]

fwrite(rows, file.path(OUTD, "table_s10_temporal_persistence.csv"))

# ---- paste-ready markdown fragment (same column set as the printed Table S10) ----
fmt_n   <- function(x) ifelse(is.na(x), "—", as.character(x))
fmt_pct <- function(x) ifelse(is.na(x), "—",
                              ifelse(x == round(x) & x %% 1 == 0, sprintf("%.0f", x), sprintf("%.1f", x)))
trans <- ifelse(is.na(rows$visit_from), paste0("(baseline) → ", rows$visit_to),
                paste0(rows$visit_from, " → ", rows$visit_to))
md <- c(
  "| Compartment | Visit transition | Sig. (preceding) | Sig. (current) | Retained | Newly sig. | Retained (%) |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %s | %s | %s | %d | %d | %d | %s |",
          rows$compartment, trans, fmt_n(rows$n_sig_preceding), rows$n_sig_current,
          rows$n_retained, rows$n_newly_sig, fmt_pct(rows$pct_retained)))
writeLines(md, file.path(OUTD, "table_s10_temporal_persistence.md"))

cat("wrote", file.path(OUTD, "table_s10_temporal_persistence.csv"), "\n")
print(rows[, .(compartment, visit_from, visit_to, n_sig_preceding, n_sig_current,
               n_retained, n_newly_sig, pct_retained)])
cat("\n--- markdown ---\n"); cat(md, sep = "\n"); cat("\n")
