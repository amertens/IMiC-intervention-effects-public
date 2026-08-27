# run-proteomics-go.R — build Table S6 (proteome GO-BP over-representation) from
# the UniProt-native result, so Table S6 and Fig 6C share ONE method.
#
# HISTORY / WHY THIS IS A THIN TRANSFORM NOW
#   This script used to run its own gene-level `enrichGO` after mapping UniProt ->
#   Entrez. That path is statistically inferior for this assay: bitr() expands one
#   measured protein into several Entrez genes and inflates the hypergeometric
#   counts, and it pooled a single global universe across studies. Fig 6C is drawn
#   from the correct UniProt-native rerun (`src/2 analysis/55-proteomics-go-uniprot.R`,
#   each protein counted once, per-cell measured background). To keep the table and
#   its own figure consistent, Table S6 is now DERIVED from that same UniProt output
#   rather than recomputed by a divergent method. The Entrez path is retired.
#
# INPUT   results/proteomics_go_uniprot.csv        (from 55-proteomics-go-uniprot.R)
# OUTPUT  results/metaboanalyst/proteomics_go/proteomics_go_tableS6.csv
#         (same filename + column schema downstream consumers already expect:
#          figure-pathway-replication-matrix.R and the manuscript Table S6.)
#
# Table S6 inclusion mirrors the submitted table's breadth (the enrichGO default:
# nominal p < 0.05 AND Storey q < 0.20) but computed with the correct protein-level
# test. The `fdr` column carries MetaboAnalyst/clusterProfiler's BH p.adjust so the
# significance frontier (fdr < 0.05) used by Fig 6C and the replication matrix is
# transparent in the table.
suppressMessages({ library(dplyr); library(readr) })

UNIPROT_CSV <- "results/proteomics_go_uniprot.csv"
OUT_DIR     <- "results/metaboanalyst/proteomics_go"
OUT_CSV     <- file.path(OUT_DIR, "proteomics_go_tableS6.csv")

P_CUTOFF <- 0.05    # nominal p (enrichGO default)
Q_CUTOFF <- 0.20    # Storey q (enrichGO default qvalueCutoff)

build_tableS6 <- function(uniprot_csv = UNIPROT_CSV, write = TRUE) {
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

if (sys.nframe() == 0) invisible(build_tableS6(write = TRUE))
