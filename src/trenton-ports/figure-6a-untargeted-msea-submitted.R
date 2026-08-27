# =============================================================================
# figure-6a-untargeted-msea-submitted.R
# Faithful R port of the SUBMITTED Fig 6A figure code, found in Trenton's older
# "imicPaperUntargetedMetabolomics.Rmd" (2025-07-26) -- the generator the 2026
# analysis notebooks lacked.
#
# SUBMITTED CONSTRUCTION (his exact pipeline):
#   For each of 14 cells (study x direction x timepoint):
#     1. significant untargeted milk features: sig==1, est sign, biomarker starts
#        "Rlc" (relabelled "rLC"), measure=="ATE"; ELICIT & Mumta-LW 1.5mo from the
#        COMBINED-arms RDS, MISAME & Mumta-LW 2mo from the STRATIFIED-arms RDS.
#     2. annotate each feature to a compound NAME via his annotation keys
#        (misame3 manual annotation + misame3_name_map_1..12 for MISAME;
#         child_elicit_mumptalw_final_annotation_key for ELICIT/Mumta-LW).
#     3. upload the distinct names to MetaboAnalyst Enrichment (ORA) with the
#        untargeted reference metabolome; download msea_ora_result_<cell>.xlsx.
#   Figure: bind all cells, enrichment_ratio = hits/expected, signed by direction
#     (down negative); grey if Raw p >= 0.05; jittered; label hits>=1 & Raw p<0.05.
#     x = signed Enrichment Ratio, y = -log10(Raw p); P<0.05 (red) & ER=2 (blue)
#     lines. (his "★Untargeted Metabolites Figure", saved for_andrew_figure_4.png.)
#
# AUTOMATION: step 3 (manual metaboanalyst.ca ORA upload/download) is replaced by
# the repo's validated run_ora() engine, so no msea_ora_result_*.xlsx download is
# needed -- the ORA is reproduced offline and fed to the exact figure code below.
#
# INPUT DEPENDENCY (why this does not run end-to-end in-repo yet): step 2 needs
# Trenton's untargeted ANNOTATION KEYS, which are NOT in the repo:
#   - misame3 manual annotation: untargeted_to_annotate_manually_annotated.xlsx
#   - misame3_name_map_1..12.csv  (MetaboAnalyst ID-conversion of the above)
#   - child_elicit_mumptalw_final_annotation_key (built from milkAnnotationKey*.csv
#     + IMiC_Azad_Metabolite_ID_Share.xlsx + child_elicit_mumptalw_three_metabolites.csv)
#   - the untargeted reference metabolome (misame3_name_map_all + the CEM name map).
# `annotate_features()` and `.untargeted_reference()` below are the hooks; supply
# those files (or point ANNOT_DIR at them) and this reproduces the submitted 6A 1:1.
# Until then, the repo's own in-repo reconstruction (fig6A-untargeted-msea.R via
# run-untargeted-msea.R, curated annotations + whole-metabolome background) is the
# working panel; note it differs from this submitted construction.
#
# Run from repo root:  Rscript src/trenton-ports/figure-6a-untargeted-msea-submitted.R
# =============================================================================
suppressMessages({ library(dplyr); library(stringr); library(ggplot2); library(ggrepel) })
source("src/metaboanalyst/R/run-ora.R")     # validated Enrichment Analysis (ORA)
source("src/metaboanalyst/R/harvest.R")

COMBINED_RDS   <- "results/combined_intervention_effects_results_combined_arms.RDS"
STRATIFIED_RDS <- "results/combined_intervention_effects_results_stratified_arms.RDS"
ANNOT_DIR      <- Sys.getenv("IMIC_UNTARGETED_ANNOT_DIR", "trenton scripts/untargeted annotation")
OUT            <- "results/trenton-ports/figure_6a"
PNG            <- "figures/trenton_ports/figure_6a_untargeted_msea_submitted.png"

# His 14 cells: (study, direction, timepoint, studytime, arm_set). ELICIT + Mumta 1.5mo
# use combined arms; MISAME + Mumta 2mo use stratified (his per-cell notes).
CELLS <- tibble::tribble(
  ~study,  ~direction,      ~timepoint,    ~studytime,             ~arms,
  "elicit","downregulated", "1 mo",        "Elicit (1 mo.)",       "combined",
  "elicit","upregulated",   "1 mo",        "Elicit (1 mo.)",       "combined",
  "elicit","downregulated", "5 mo",        "Elicit (5 mo.)",       "combined",
  "elicit","upregulated",   "5 mo",        "Elicit (5 mo.)",       "combined",
  "misame","downregulated", "14-21 days",  "Misame (14-21 days)",  "stratified",
  "misame","upregulated",   "14-21 days",  "Misame (14-21 days)",  "stratified",
  "misame","downregulated", "1-2 months",  "Misame (1-2 mo.)",     "stratified",
  "misame","upregulated",   "1-2 months",  "Misame (1-2 mo.)",     "stratified",
  "misame","downregulated", "3-4 months",  "Misame (3-4 mo.)",     "stratified",
  "misame","upregulated",   "3-4 months",  "Misame (3-4 mo.)",     "stratified",
  "vital", "downregulated", "1.5 mo",      "Vital (1.5 mo.)",      "combined",
  "vital", "upregulated",   "1.5 mo",      "Vital (1.5 mo.)",      "combined",
  "vital", "downregulated", "2 mo",        "Vital (2 mo.)",        "stratified",
  "vital", "upregulated",   "2 mo",        "Vital (2 mo.)",        "stratified")

# Reproducible: significant untargeted features for a cell (before annotation).
sig_untargeted <- function(dat, studytime_value, direction) {
  dat %>% filter(studytime == studytime_value, contrast != "control",
                 str_starts(biomarker, "Rlc"), sig == 1, measure == "ATE",
                 if (direction == "upregulated") est > 0 else est < 0) %>%
    mutate(biomarker = str_replace(biomarker, "^Rlc", "rLC"))
}

# HOOK: map feature ids (biomarker) -> compound names using Trenton's annotation
# keys. Requires the files listed in the header; hard-fails clearly if absent.
annotate_features <- function(sig_df, study) {
  stop("figure-6a port: untargeted annotation keys not found in '", ANNOT_DIR, "'. ",
       "Supply Trenton's untargeted annotation files (see header) to reproduce the ",
       "submitted Fig 6A end-to-end; until then use the in-repo reconstruction at ",
       "figure-scripts/manuscript_figures/fig6A-untargeted-msea.R.", call. = FALSE)
}
.untargeted_reference <- function() {
  stop("figure-6a port: untargeted reference metabolome not available (see header).", call. = FALSE)
}

# --- the SUBMITTED figure code (his ★Untargeted Metabolites Figure), verbatim in
#     logic. Input: one row per (cell x pathway) with pathway,total,expected,hits,
#     `Raw p`,direction,study,timepoint. -------------------------------------------
plot_untargeted_msea <- function(msea_all, out_png = PNG) {
  tableau20 <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1",
                 "#FF9DA7","#9C755F","#BAB0AC","#86BCD6","#FFBE7D","#FF5850","#A0CBE8",
                 "#8CD17D","#B6992D","#499894","#FABFD2","#D37295","#B7B7B7")
  set.seed(123)
  d <- msea_all %>%
    mutate(study = recode(study, elicit = "ELICIT", misame = "MISAME-III", vital = "MumptaLW", .default = study),
           enrichment_ratio  = hits / expected,
           enrichment_signed = if_else(tolower(direction) == "downregulated", -enrichment_ratio, enrichment_ratio),
           is_grey = (-log10(`Raw p`)) < -log10(0.05),
           study_color = ifelse(is_grey, "Not Significant", study),
           jitter_x = ifelse(is_grey, enrichment_signed, enrichment_signed + runif(n(), -1, 1)),
           jitter_y = ifelse(is_grey, -log10(`Raw p`), pmax(-log10(0.05), -log10(`Raw p`) + runif(n(), -0.3, 0.3))))
  studies <- unique(d$study); cols <- setNames(tableau20[seq_along(studies)], studies)
  cols <- c("Not Significant" = "grey", cols)
  p <- ggplot(d, aes(jitter_x, jitter_y)) +
    geom_point(aes(color = study_color), size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 2, linetype = "dashed", color = "blue") +
    geom_label_repel(data = subset(d, hits >= 1 & `Raw p` < 0.05), aes(label = pathway, color = study_color),
                     size = 2, show.legend = FALSE, max.overlaps = 200) +
    scale_color_manual(values = cols) +
    labs(x = "Enrichment Ratio (signed)", y = expression(-Log[10]*"(Raw P)"), color = "Study") +
    scale_x_continuous(limits = c(-10, 30), breaks = seq(-10, 30, 5)) +
    scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 0.5)) +
    theme_classic(base_size = 7) + theme(legend.position = "bottom")
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_png, p, width = 210/2, height = 297/3, units = "mm", dpi = 600, device = ragg::agg_png)
  cat("wrote", out_png, "\n"); invisible(p)
}

# Reproduce each cell's ORA via run_ora on the annotated significant names, then plot.
if (sys.nframe() == 0) {
  combined <- readRDS(COMBINED_RDS); stratified <- readRDS(STRATIFIED_RDS)
  ref <- .untargeted_reference()                       # hard-fails until keys supplied
  rows <- list()
  for (i in seq_len(nrow(CELLS))) {
    ce  <- CELLS[i, ]
    dat <- if (ce$arms == "combined") combined else stratified
    names_vec <- sig_untargeted(dat, ce$studytime, ce$direction) %>% annotate_features(ce$study)
    if (!length(names_vec)) next
    mSet <- tryCatch(run_ora(names_vec, reference_names = ref), error = function(e) NULL)
    if (is.null(mSet)) next
    res <- harvest_results(mSet) %>% transmute(pathway, total, expected, hits, `Raw p` = raw_p,
                                               direction = ce$direction, study = ce$study, timepoint = ce$timepoint)
    rows[[length(rows) + 1]] <- res
  }
  msea_all <- bind_rows(rows)
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(msea_all, file.path(OUT, "figure_6a_untargeted_msea_all_cells.csv"))
  plot_untargeted_msea(msea_all)
}
