# label-map.R, curated compound-name synonyms -> MetaboAnalyst-recognized names.
# Source of truth, extracted verbatim from "trenton scripts/2. Scripts/Primary Outcomes.Rmd".
suppressMessages(library(stringr))

LABEL_MAP <- c(
  "nicotinamide mononucleotide"                   = "Nicotinic acid mononucleotide",
  "vitamin B3 \\( expressed as NAM\\)"            = "Niacinamide",
  "nudifloramide \\(niacin catabolite\\)"         = "Nudifloramide",
  "free thiamin"                                  = "Thiamine",
  "total vitamin B1 \\(expressed as thiamin\\)"   = "Thiamine",
  "vitamin B2 \\(expressed as riboflavin\\)"      = "Riboflavin",
  "vitamin B12 \\(expressed as cyanocobalamin\\)" = "Cyanocobalamin",
  "vitamin B6 \\(expressed as PL\\)"              = "Pyridoxine",
  "^a\\.tocopherol$|^A\\.tocopherol$"             = "Alpha-tocopherol",
  "^g\\.tocopherol$|^G\\.tocopherol$"             = "Gamma-tocopherol",
  "^vitamin\\.a$|^Vitamin\\.a$"                   = "Vitamin A"
)

# Apply the curated rename map to a vector of raw feature labels (label_f).
apply_label_map <- function(label_f) {
  str_replace_all(label_f, LABEL_MAP)
}
