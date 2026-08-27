# =============================================================================
# blood_compartment_volcano_plots.R
#
# Volcano plots for the blood-compartment intervention-effect results
# (maternal plasma / VAMS prenatal / VAMS postnatal maternal+infant / proteomics
# depleted+naive), for both the arm-stratified and combined-arms codings.
#
# x = ATE, y = -log10(raw p); points coloured by significance (raw vs FDR), using
# the repo's plot_imic_volcano(). One figure per dataset, faceted visit x contrast.
#
# Reads the cleaned results written by clean_blood_results.R. Tag-aware:
# set BLOOD_ADJUST=TRUE to plot the adjusted results instead of unadjusted.
# Run AFTER the analysis + clean step. Output: figures/blood_volcano/
# =============================================================================

if (!exists("BLOOD_ORCHESTRATED")) { rm(list = ls()); source(paste0(here::here(), "/src/0-config.R")) }
if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- FALSE
.tag <- if (BLOOD_ADJUST) "adjusted_" else ""

outdir <- paste0(here::here(), "/figures/blood_volcano")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

clean_path <- function(kind)  # kind = "" (stratified) or "combined_arms_"
  paste0(here::here(), "/results/blood_compartment_", .tag, kind,
         "intervention_effects_results_clean.RDS")

# plot_imic_volcano() needs pval/pval_adj/est + a label_f column; label FDR-sig
# features by their biomarker id (metabolite m/z-RT or UniProt protein).
prep_volcano <- function(path) {
  if (!file.exists(path)) { warning("missing: ", path); return(NULL) }
  readRDS(path) %>%
    dplyr::ungroup() %>%
    filter(measure == "ATE") %>%
    mutate(label_f = biomarker)
}

# One faceted volcano (visit x contrast) per dataset; size adapts to the grid.
save_volcanos <- function(d, coding) {
  if (is.null(d) || nrow(d) == 0) return(invisible())
  for (ds in unique(d$dataset)) {
    di <- d %>% filter(dataset == ds)
    p <- plot_imic_volcano(di, title = paste0(ds, "  (", coding, " arms, ",
                                              if (.tag == "") "unadjusted" else "adjusted", ")"),
                           facet_type = "arm", labels = TRUE, label_type = "label",
                           overlap_n = 15)
    w <- 3 + 3.2 * max(1, dplyr::n_distinct(di$contrast))
    h <- 2 + 2.2 * max(1, dplyr::n_distinct(di$visit))
    ggsave(p, file = paste0(outdir, "/blood_volcano_", .tag, coding, "_", ds, ".png"),
           width = w, height = h, limitsize = FALSE)
  }
  message("[", coding, "] wrote ", dplyr::n_distinct(d$dataset), " volcano figures to ", outdir)
}

save_volcanos(prep_volcano(clean_path("")),               "stratified")
save_volcanos(prep_volcano(clean_path("combined_arms_")), "combined")
