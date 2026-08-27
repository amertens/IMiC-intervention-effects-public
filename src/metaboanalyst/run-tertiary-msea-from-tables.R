# run-tertiary-msea-from-tables.R -- FAITHFUL reproduction of the SUBMITTED Fig 5B.
#
# The submitted tertiary MSEA panel was drawn from Trenton's downloaded
# MetaboAnalyst ORA result tables (data/msea/msea_ora_result_significant_*_tert_*_
# hmdb*.csv), exactly as Fig 6A is reproduced from his untargeted exports via
# run-untargeted-msea.R. This script is the tertiary analogue: it reads those
# per-cell tables directly (no re-run), so the panel matches the SUBMISSION
# (including its wide enrichment-ratio scale, which comes from his background N).
#
# NOTE: this is the "closest repo replication of the exact submitted 5B". The
# reproducible-pipeline alternative (re-run ORA via MetaboAnalystR with a named
# reference metabolome) is run-tertiary-msea.R -- same pathways, cleaner method,
# but different points/scale because the query construction and background differ.
#
# Each table: rows = SMPDB metabolite sets; cols = total, expected, hits, Raw p,
# Holm p, FDR. We reshape to the panel schema and sign the enrichment ratio by
# direction (down = negative). Metabolite tables only (lipid-set tables excluded);
# where a cell has both a bare `_hmdb` and a `_hmdb_metabolites` file, the
# `_metabolites` file wins.
#
# Run from repo root:  Rscript src/metaboanalyst/run-tertiary-msea-from-tables.R
suppressMessages({ library(dplyr); library(readr); library(stringr); library(tibble) })

MSEA_DIR <- "data/msea"
OUT      <- "results/metaboanalyst/tertiary_msea"

.STUDY    <- c(elicit = "ELICIT", misame = "MISAME-III", vital = "Mumta-LW")
.CONTRAST <- c(elicit = "Nico",   misame = "BEP",        vital = "BEP")
.TP       <- c("1mo" = "1 mo.", "5mo" = "5 mo.", "1421" = "14-21 days",
               "12mo" = "1-2 mo.", "34mo" = "3-4 mo.",
               "15mo" = "1.5 mo.", "2mo" = "2 mo.")

parse_cell <- function(fn) {
  b <- sub("\\.csv$", "", sub("^msea_ora_result_significant_", "", fn))
  study_key <- str_extract(b, "^(elicit|misame|vital)")
  b2  <- sub("^(elicit|misame|vital)_tert_", "", b)
  dir_key <- str_extract(b2, "^(down|up)regulated")
  # timepoint token = the piece before _combined/_stratified
  tp_key  <- str_extract(sub("^(down|up)regulated_", "", b2), "^[0-9a-z]+")
  list(study     = unname(.STUDY[study_key]),
       contrast  = unname(.CONTRAST[study_key]),
       direction = if (identical(dir_key, "downregulated")) "down" else "up",
       timepoint = unname(.TP[tp_key]),
       cell_base = sub("_hmdb(_metabolites|_lipids)?$", "", b))  # for metabolites-vs-bare dedup
}

run_tertiary_msea_from_tables <- function(write = TRUE) {
  files <- list.files(MSEA_DIR, pattern = "_tert_.*_hmdb.*\\.csv$")
  # Keep only well-formed names (drops a doubly-prefixed elicit lipids export that
  # otherwise parses to study/timepoint = NA).
  files <- files[grepl("^msea_ora_result_significant_(elicit|misame|vital)_tert_", files)]
  meta  <- lapply(files, parse_cell)
  # The submitted 5B combined BOTH the HMDB metabolite-set AND lipid-set enrichment
  # runs: fatty-acid pathways (Oxidation of BCFAs, Carnitine Synthesis, Ketone Body)
  # live in the *_lipids tables and carry the large negative enrichment ratios (< -250)
  # seen in the submitted panel. Keep every *_metabolites and *_lipids table; use a bare
  # *_hmdb table only for a cell that has no *_metabolites table (e.g. elicit 5 mo).
  base   <- vapply(meta, `[[`, "", "cell_base")
  is_met <- grepl("_hmdb_metabolites\\.csv$", files)
  is_lip <- grepl("_hmdb_lipids\\.csv$", files)
  keep   <- vapply(seq_along(files), function(i)
    is_met[i] || is_lip[i] || !any(is_met & base == base[i]), logical(1))
  files <- files[keep]; meta <- meta[keep]
  message("reading ", length(files), " tertiary ORA tables from ", MSEA_DIR)

  rows <- Map(function(fn, ce) {
    x <- read.csv(file.path(MSEA_DIR, fn), row.names = 1, check.names = FALSE,
                  stringsAsFactors = FALSE)
    if (!nrow(x)) return(NULL)
    num <- function(col) suppressWarnings(as.numeric(x[[col]]))
    tibble(study = ce$study, timepoint = ce$timepoint, contrast = ce$contrast,
           direction = ce$direction, pathway = rownames(x),
           total = num("total"), expected = num("expected"), hits = num("hits"),
           raw_p = num("Raw p"), fdr_native = num("FDR"),
           enrichment_ratio = ifelse(ce$direction == "down", -1, 1) * (num("hits") / num("expected")))
  }, files, meta)
  tab <- bind_rows(rows) %>% arrange(raw_p)

  if (isTRUE(write)) {
    dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
    write_csv(tab, file.path(OUT, "tertiary_msea_fromTables.csv"))
    write_csv(filter(tab, fdr_native < 0.05),
              file.path(OUT, "tertiary_msea_fromTables_perContrastFDRsig_tableS4.csv"))
  }
  tab
}

if (sys.nframe() == 0) {
  tab <- run_tertiary_msea_from_tables(write = TRUE)
  message("done tertiary MSEA (from tables): ", nrow(tab), " pathway rows across ",
          dplyr::n_distinct(tab[c("study","timepoint","direction")]), " cells; ",
          sum(tab$fdr_native < 0.05, na.rm = TRUE), " FDR-significant; ",
          sum(tab$raw_p < 0.05, na.rm = TRUE), " nominally significant")
}
