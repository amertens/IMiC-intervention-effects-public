# =============================================================================
# 41-fat-synthesis-timepoint-table.R
#
# Tabulates the signed direction of each milk lipid / fat-synthesis Mummichog
# pathway at all three MISAME-III milk visits (14-21 d, 1-2 mo, 3-4 mo) side by side.
# The cross-compartment comparison uses the 1-2 mo milk visit (when all sample types
# were collected), so this checks that a direction seen there also holds at the
# first and third visits rather than reflecting one timepoint (for example, milk
# volume dilution at 1-2 mo). It backs the Discussion statement that de novo
# fatty-acid biosynthesis was decreased in milk. The direction measure is the
# per-visit binomial sign test from script 33 (Mummichog enrichment alone is
# direction-agnostic); this script only pivots that long table to one row per
# pathway with timepoint columns.
#
# Input  : results/signed_pathway_direction_milk.csv (script 33)
# Output : results/fat_synthesis_timepoint_table.csv
# [needs restricted data] (upstream inputs of script 33)
# =============================================================================
suppressMessages({library(data.table)})
root <- paste0(here::here(), "/")

inf <- paste0(root, "results/signed_pathway_direction_milk.csv")
if (!file.exists(inf)) stop("Run 33-signed-pathway-direction-milk.R first to create: ", inf)
d <- fread(inf)

# lipid / fat-synthesis pathways (de novo synthesis is the pathway of interest; the
# related lipid pathways show whether the direction is coherent).
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
cat("\nDe novo fatty acid biosynthesis:\n")
print(fs[grepl("de novo fatty acid", pathway, ignore.case = TRUE),
         .(visit, n_members, frac_up, mean_signed, direction, sign_test_p)][order(visit)])
cat("\nRead: a pathway 'down' at 1-2 mo that is also 'down' at 14-21 d and 3-4 mo is\n",
    "robust to the single-timepoint matching concern; a direction that appears only\n",
    "at 1-2 mo should be interpreted with the milk-volume-dilution caveat.\n")
cat("\nSaved: results/fat_synthesis_timepoint_table.csv\n")
