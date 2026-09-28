# =============================================================================
# run-proteomics-go.R
#
# Builds Table S7 (Gene Ontology biological-process over-representation of the
# untargeted milk proteome) by filtering the UniProt-level enrichment written by
# src/2 analysis/55-proteomics-go-uniprot.R (clusterProfiler::enricher with a
# UniProt-to-GO map; each protein counted once, per-cell measured background).
# Fig 6C reads that same file, so the table and figure share one method.
# Inclusion follows the enrichGO defaults (nominal P < 0.05 and Storey q < 0.20);
# the `fdr` column is the BH-adjusted P used as the significance threshold
# (fdr < 0.05) in Fig 6C.
#
# Inputs : results/proteomics_go_uniprot.csv  (src/2 analysis/55-proteomics-go-uniprot.R)
# Outputs: results/metaboanalyst/proteomics_go/proteomics_go_pathways.csv  (Table S7)
#
# Run from repo root: Rscript src/metaboanalyst/run-proteomics-go.R
# =============================================================================
suppressMessages({ library(dplyr); library(readr) })

UNIPROT_CSV <- "results/proteomics_go_uniprot.csv"
OUT_DIR     <- "results/metaboanalyst/proteomics_go"
OUT_CSV     <- file.path(OUT_DIR, "proteomics_go_pathways.csv")

P_CUTOFF <- 0.05    # nominal p (enrichGO default)
Q_CUTOFF <- 0.20    # Storey q (enrichGO default qvalueCutoff)

build_go_table <- function(uniprot_csv = UNIPROT_CSV, write = TRUE) {
  if (!file.exists(uniprot_csv)) {
    stop("UniProt-native proteomics result not found at '", uniprot_csv, "'. ",
         "Run src/2 analysis/55-proteomics-go-uniprot.R first (it writes this file).",
         call. = FALSE)
  }
  u <- read_csv(uniprot_csv, show_col_types = FALSE)

  # qvalue can be absent/all-NA in edge cases; fall back to BH p.adjust < 0.05 so
  # the table is never silently emptied by a missing Storey q.
  has_q <- "qvalue" %in% names(u) && any(!is.na(u$qvalue))
  keep <- if (has_q) (u$pvalue < P_CUTOFF & u$qvalue < Q_CUTOFF)
          else       (u$p.adjust < P_CUTOFF)
  keep[is.na(keep)] <- FALSE

  tab <- u[keep, , drop = FALSE] %>%
    transmute(
      ID,
      Description,
      gene_ratio      = GeneRatio,
      fold_enrichment = signed_fold_enrichment,   # signed by direction (down = negative)
      p_value         = pvalue,
      fdr             = p.adjust,                  # BH-adjusted (the significance metric)
      study,
      timepoint,
      contrast,
      regulation      = direction                 # "Upregulated" / "Downregulated"
    ) %>%
    arrange(p_value)

  if (isTRUE(write)) {
    dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
    write_csv(tab, OUT_CSV)
    message("wrote ", OUT_CSV, "  (", nrow(tab), " GO-BP terms; ",
            sum(tab$fdr < 0.05, na.rm = TRUE), " at FDR<0.05)")
  }
  tab
}

if (sys.nframe() == 0) invisible(build_go_table(write = TRUE))
