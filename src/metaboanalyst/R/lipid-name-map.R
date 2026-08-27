# lipid-name-map.R -- convert Quant 500 verbose lipid nomenclature to the LIPID
# MAPS abbreviations that MetaboAnalystR's lipid_compound_db (CrossReferencing
# with lipid = TRUE) recognises.
#
# WHY: MetaboAnalystR does NOT run the MetaboAnalyst *web* smart-matcher that
# normalises vendor lipid names, so our verbose Quant 500 names match only ~2%
# offline. Reformatted to LIPID MAPS abbreviations they match ~57% (measured on
# the 469 tertiary lipid names). The residual ~43% are Quant 500 triglycerides
# named as "one acyl chain + remainder" (e.g. "Triacylglyceride (18:1_30:2)",
# where 30:2 is the summed remainder of the other two chains) -- that hybrid is
# not a standard LIPID MAPS species and cannot resolve to a compound record.
#
# This is the offline stand-in for the web "feature type = lipids" step. Feed the
# converted names to run_ora(..., lipid = TRUE).
suppressMessages(library(stringr))

# Tertiary Quant 500 categories that are lipids (-> the lipid DB name-matching).
LIPID_CATEGORIES <- c(
  "Triglycerides", "Diglycerides", "Ceramides", "Hexosylceramides",
  "Dihexosylceramides", "Trihexosylceramides", "Dihydroceramides",
  "Phosphatidylcholines", "Lysophosphatidylcholines", "Cholesteryl Esters",
  "Sphingomyelins")
  # NB: acylcarnitines / free fatty acids are kept in the METABOLITE pass. Moving
  # them to the lipid feature-type run (to match the submitted's split) was tested
  # and made recovery WORSE (16 -> 15 of 20; it dropped Oxidation of BCFAs), so the
  # 4 borderline pathways are not recoverable by a feature-type split.

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
