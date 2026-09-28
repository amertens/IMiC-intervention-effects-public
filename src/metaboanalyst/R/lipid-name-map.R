# =============================================================================
# lipid-name-map.R
#
# Converts Quant 500 lipid names to the LIPID MAPS abbreviations that
# MetaboAnalystR's lipid_compound_db recognises (CrossReferencing with
# lipid = TRUE); the scripted stand-in for the web tool's "feature type = lipids"
# step. Used by run-tertiary-msea-dual.R's lipid pass (Fig 5B, Table S3).
# MetaboAnalystR does not include the web tool's lipid-name normaliser, so the
# verbose Quant 500 names match only about 2% of the time; the LIPID MAPS forms
# match about 57% (469 tertiary lipid names). Most of the rest are triglycerides
# named as "one acyl chain + remainder" (e.g. "Triacylglyceride (18:1_30:2)"),
# which is not a standard LIPID MAPS species and has no compound record.
#
# Inputs : none
# Outputs: none (defines LIPID_CATEGORIES and to_lipidmaps())
# =============================================================================
suppressMessages(library(stringr))

# Tertiary Quant 500 categories that are lipids (-> the lipid DB name-matching).
# Acylcarnitines and free fatty acids stay in the metabolite pass: moving them to
# the lipid pass recovered fewer of the submitted Fig 5B pathways (15 instead of
# 16 of 20; it lost oxidation of branched-chain fatty acids).
LIPID_CATEGORIES <- c(
  "Triglycerides", "Diglycerides", "Ceramides", "Hexosylceramides",
  "Dihexosylceramides", "Trihexosylceramides", "Dihydroceramides",
  "Phosphatidylcholines", "Lysophosphatidylcholines", "Cholesteryl Esters",
  "Sphingomyelins")

to_lipidmaps <- function(x) {
  x <- str_replace(x, "^Triacylglyceride \\((.*)\\)$",  "TG(\\1)")
  x <- str_replace(x, "^Diacylglyceride \\((.*)\\)$",   "DG(\\1)")
  x <- str_replace(x, "^Dihydroceramide \\((.*)\\)$",   "Cer(\\1)")   # d18:0 backbone
  x <- str_replace(x, "^Trihexosylceramide \\((.*)\\)$","Hex3Cer(\\1)")
  x <- str_replace(x, "^Dihexosylceramide \\((.*)\\)$", "Hex2Cer(\\1)")
  x <- str_replace(x, "^Hexosylceramide \\((.*)\\)$",   "HexCer(\\1)")
  x <- str_replace(x, "^Ceramide \\((.*)\\)$",          "Cer(\\1)")
  x <- str_replace(x, "^Cholesteryl ester (\\d+:\\d+)$","CE(\\1)")
  x <- str_replace(x, "^Phosphatidylcholine a[ae] C(\\d+:\\d+)$", "PC(\\1)")
  x <- str_replace(x, "^Lysophosphatidylcholine a C(\\d+:\\d+)$", "LPC(\\1)")
  x <- str_replace(x, "^Hydroxysphingomyelin C(\\d+:\\d+)$", "SM(\\1)")
  x <- str_replace(x, "^Sphingomyelin C(\\d+:\\d+)$",   "SM(\\1)")
  x
}
