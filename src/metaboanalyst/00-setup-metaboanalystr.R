# =============================================================================
# 00-setup-metaboanalystr.R
#
# One-time environment setup for the scripts in src/metaboanalyst/: installs the
# CRAN and Bioconductor dependencies, then MetaboAnalystR from GitHub
# (xia-lab/MetaboAnalystR). The published results were produced with
# MetaboAnalystR 4.3.0 under R 4.4.2; this installs the current GitHub version, so
# compare against the recorded session information if results differ. On Windows,
# MetaboAnalystR needs Rtools44 to compile. Safe to re-run: every install step is
# guarded by requireNamespace().
#
# Inputs : none
# Outputs: src/metaboanalyst/env/sessionInfo.txt  (R and package versions)
#
# Run from repo root: Rscript src/metaboanalyst/00-setup-metaboanalystr.R
# =============================================================================

message("R version: ", getRversion())
if (.Platform$OS.type == "windows" && !pkgbuild::has_rtools(debug = FALSE)) {
  warning("Rtools not detected. MetaboAnalystR needs Rtools44 to compile. ",
          "Install from https://cran.r-project.org/bin/windows/Rtools/ then re-run.")
}

# --- CRAN helpers -----------------------------------------------------------
# Install each only if it is not already available (requireNamespace guard).
cran_pkgs <- c("devtools", "pkgbuild", "remotes", "qs", "crmn", "httr",
               "jsonlite", "tidyverse", "testthat")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
}

# --- Bioconductor dependency chain -----------------------------------------
# MetaboAnalystR depends on many Bioconductor packages that are not on CRAN;
# install any that are missing in a single BiocManager call.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz",
               "preprocessCore", "siggenes", "multtest", "RBGL", "edgeR",
               "fgsea", "sva", "limma", "KEGGgraph", "SSPA")
missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

# --- MetaboAnalystR itself ---------------------------------------------------
if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
  devtools::install_github("xia-lab/MetaboAnalystR",
                           build = TRUE, build_vignettes = FALSE, build_manual = FALSE)
}

# --- Record environment for reproducibility ---------------------------------
dir.create("src/metaboanalyst/env", showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), "src/metaboanalyst/env/sessionInfo.txt")
message("MetaboAnalystR installed: ",
        requireNamespace("MetaboAnalystR", quietly = TRUE))
