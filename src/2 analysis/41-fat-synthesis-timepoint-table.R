# =============================================================================
# 41-fat-synthesis-timepoint-table.R
#
# Action item (2026-06-30 IMiC meeting): "add first and third timepoint columns to
# the fat-synthesis analysis." The cross-compartment work was matched at the second
# milk timepoint (1-2 mo, when all sample types were taken), and Andrew flagged the
# risk of over-interpreting a single-timepoint fat-synthesis result (possible milk-
# volume dilution at 1-2 mo). This builds the robustness table: for each lipid /
# fat-synthesis milk pathway, the signed direction at ALL THREE milk timepoints
# (14-21 d, 1-2 mo, 3-4 mo) side by side, so the reader can see whether the
# second-timepoint direction holds at the first and third.
#
# Source = the per-visit signed sign-test already produced by
# 33-signed-pathway-direction-milk.R (the defensible DIRECTIONAL measure -- mummichog
# enrichment alone is direction-agnostic, see response-to-reviewers 2.4). This script
# just pivots that long table to one row per pathway with timepoint columns.
#
# Out: results/fat_synthesis_timepoint_table.csv
# =============================================================================
suppressMessages({library(data.table)})
root <- paste0(here::here(), "/")

inf <- paste0(root, "results/signed_pathway_direction_milk.csv")
if (!file.exists(inf)) stop("Run 33-signed-pathway-direction-milk.R first to create: ", inf)
d <- fread(inf)

# lipid / fat-synthesis pathways (de novo synthesis is the headline; the related
# lipid pathways are included so the table shows whether the direction is coherent).
# Keyword OR-match on the pathway name; require >=4 members so a direction is meaningful.
LIPID <- paste0("fatty acid|arachidon|leukotriene|carnitine|linole|omega|",
                "glycerophospho|unsaturated|prostaglandin|sphingolipid|lipid")
fs <- d[grepl(LIPID, pathway, ignore.case = TRUE) & n_members >= 4]

# enforce timepoint order T1 -> T2 -> T3 (factor levels control the wide column order)
vlev <- c("14-21 days", "1-2 mo.", "3-4 mo.")
fs[, visit := factor(visit, levels = vlev)]
fs <- fs[!is.na(visit)]
if (!nrow(fs)) stop("No lipid pathways with >=4 members found across milk timepoints.")

# pivot long -> wide: one row per pathway, one column block per timepoint, so the
# three timepoints sit side by side for each reported quantity.
wide <- dcast(fs, pathway ~ visit,
              value.var = c("direction", "frac_up", "mean_signed", "sign_test_p", "n_members"))
setorder(wide, pathway)
fwrite(wide, paste0(root, "results/fat_synthesis_timepoint_table.csv"))

# console: compact direction matrix + the de novo FA biosynthesis detail
cat("=== Fat-synthesis / lipid milk pathways: signed direction across timepoints ===\n")
print(dcast(fs, pathway ~ visit, value.var = "direction"))
cat("\nDe novo fatty acid biosynthesis (the headline fat-synthesis pathway):\n")
print(fs[grepl("de novo fatty acid", pathway, ignore.case = TRUE),
         .(visit, n_members, frac_up, mean_signed, direction, sign_test_p)][order(visit)])
cat("\nRead: a pathway 'down' at 1-2 mo that is also 'down' at 14-21 d and 3-4 mo is\n",
    "robust to the single-timepoint matching concern; a direction that appears only\n",
    "at 1-2 mo should be interpreted with the milk-volume-dilution caveat.\n")
cat("\nSaved: results/fat_synthesis_timepoint_table.csv\n")
