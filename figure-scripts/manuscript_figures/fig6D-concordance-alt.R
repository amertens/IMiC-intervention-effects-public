# fig6D-concordance-alt.R
# =============================================================================
# ALTERNATE Panel D design (PRESERVED, not the current manuscript panel).
# This is the concordance dot-plot from the Andrew x Trenton call (2026-08):
# BOTH directions, one dot per compartment with 95% CI, ordered by number of
# compartments reached. The CURRENT manuscript Panel D is Trenton's
# pathway-annotated UPREGULATED slope plot -- see fig6D-crosscompartment.R
# (which ports his "Compartment Tracking.Rmd" figure web-independently). Kept here
# so the concordance framing (which also shows downregulated transfers) is not lost.
#
# cross-compartment transfer of BEP-affected metabolite features.
#
# Design (per the Andrew x Trenton call, 2026-08): each row is ONE metabolite,
# plotted at its scaled ATE with a 95% CI and a slight vertical dodge per
# compartment. ONLY FDR-significant appearances are drawn -- Trenton was explicit
# that non-significant (open) dots should NOT appear -- so a compound shows one
# filled dot per compartment where it is significant (2 to 4: Maternal plasma,
# Maternal VAMS, Milk, Infant VAMS), colored by compartment.
#
# Inclusion: metabolites FDR-significant in >= 2 compartments in the SAME
# direction (concordant transfer). Same-name adducts/features are collapsed to
# one row per compound (Option A), keeping the strongest (min-fdr) appearance per
# compartment.
#
# Labels prefer the Sapient-curated name (results/sapient_annotation_handoff.csv,
# by m/z), else Trenton's putative_name, else m/z. Non-informative placeholder
# names (e.g. "Acid") fall back to m/z.
#
# Inputs (in-repo):
#   results/compartment_tracking/compartment_named_features_slim.RDS
#   results/sapient_annotation_handoff.csv
# Output: figures/figure6_panelD_crosscompartment.png
#
# Run from repo root: Rscript figure-scripts/manuscript_figures/fig6D-crosscompartment.R
# =============================================================================
suppressMessages({ library(dplyr); library(stringr); library(ggplot2); library(tidyr) })

# Trenton's cross-compartment linkage (exact feature-id / identical normalized
# putative name / same ion mode & mass within 25 ppm), ported faithfully in
# src/trenton-ports/compartment-tracking.R. We consume ct_tracking_components()
# from it so Panel D groups features into "the same metabolite" exactly as his
# Compartment Tracking.Rmd does -- replacing the earlier m/z-to-3-decimals key,
# which was ~15x tighter than his 25 ppm window and split concordant features
# (e.g. the milk 286.198 vs infant 286.202 carnitine, ~14 ppm apart) so nothing
# reached 4 compartments. (Sourced BEFORE our own constants so ours win.)
source("src/trenton-ports/compartment-tracking.R")

SLIM <- "results/compartment_tracking/compartment_named_features_slim.RDS"
SAP  <- "results/sapient_annotation_handoff.csv"
# Milk untargeted ATE results -- the ONLY compartment whose tscore is absent from
# the slim RDS; used to recover milk 95% CIs (carries est/cil/ciu per feature).
MILK_ATE <- "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS"
OUT  <- "figures/figure6_panelD_concordance.png"   # alternate output; does NOT feed the composite

COMP_LEVELS <- c("Maternal plasma", "Maternal VAMS", "Milk", "Infant VAMS")
COMP_COLS   <- c("Maternal plasma"="#4E79A7", "Maternal VAMS"="#76B7B2",
                 "Milk"="#F28E2B", "Infant VAMS"="#E15759")
# Placeholder / non-informative annotation strings that should fall back to m/z.
GENERIC_NAMES <- c("Acid")

# Robustly coerce a `significant` column to logical regardless of how the slim RDS
# stored it (logical, "TRUE"/"FALSE" strings, 1/0 numeric, or a factor of any of
# these). The old `significant %in% c(TRUE, "TRUE")` silently returned all-FALSE
# for numeric/factor storage, which would empty the panel with no error.
.as_sig <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x == 1)
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "significant")
}

build_panelD <- function(slim_path = SLIM, sap_path = SAP, out_png = OUT) {
  # --- sapient-curated name lookup by rounded m/z -----------------------------
  sap <- utils::read.csv(sap_path, stringsAsFactors = FALSE, check.names = FALSE)
  sap$sap_name <- ifelse(!is.na(sap$SAPIENT_identity) & sap$SAPIENT_identity != "",
                         sap$SAPIENT_identity, sap$our_putative_name)
  sap_lk <- sap %>% filter(!is.na(sap_name), sap_name != "") %>%
    mutate(mz2 = sprintf("%.2f", as.numeric(mz))) %>%
    distinct(mz2, .keep_all = TRUE) %>% select(mz2, sap_name)

  d <- readRDS(slim_path) %>%
    filter(!is.na(mz), !is.na(effect_size), compartment %in% COMP_LEVELS) %>%
    mutate(mz2 = sprintf("%.2f", mz), is_sig = .as_sig(significant))

  # Guard: if the `significant` column type ever changes so that nothing is flagged,
  # the >=2-compartment concordance filter would silently yield an empty panel.
  if (!any(d$is_sig, na.rm = TRUE))
    stop("figure6-panelD: no rows flagged significant after coercion -- check the ",
         "`significant` column in ", slim_path, call. = FALSE)

  # --- Trenton's linkage: group SIGNIFICANT features into concordant components --
  # ct_tracking_components() adds `tracking_group` (connected components of his
  # exact-id / same-name / 25 ppm same-direction matches). A "compound" is a group.
  sig <- d %>% filter(is_sig) %>% ct_tracking_components() %>%
    mutate(se  = ifelse(!is.na(tscore) & tscore != 0, abs(effect_size / tscore), NA_real_),
           cil = effect_size - 1.96 * se, ciu = effect_size + 1.96 * se)

  # --- recover MILK 95% CIs ---------------------------------------------------
  # The slim RDS carries NO tscore for Milk (se -> NA -> milk dots drew with no
  # error bar). The milk untargeted ATE results DO carry cil/ciu; join them on the
  # feature id -- case-insensitive, since the slim ids are "rLC_*" and the source
  # ids are "Rlc_*" -- disambiguated by the effect size (a feature recurs across
  # visits), then overwrite the milk cil/ciu.
  milk_ci <- readRDS(MILK_ATE)
  milk_ci <- milk_ci[milk_ci$contrast == "BEP", c("biomarker", "est", "cil", "ciu")]
  milk_ci <- milk_ci %>%
    transmute(mkey = tolower(biomarker), er = round(est, 5),
              m_cil = cil, m_ciu = ciu) %>%
    distinct(mkey, er, .keep_all = TRUE)
  sig <- sig %>%
    mutate(mkey = tolower(feature), er = round(effect_size, 5)) %>%
    left_join(milk_ci, by = c("mkey", "er")) %>%
    mutate(cil = ifelse(compartment == "Milk" & !is.na(m_cil), m_cil, cil),
           ciu = ifelse(compartment == "Milk" & !is.na(m_ciu), m_ciu, ciu)) %>%
    select(-mkey, -er, -m_cil, -m_ciu)

  # concordant transfer: tracking groups significant in >= 2 distinct compartments
  incl <- sig %>% distinct(tracking_group, compartment) %>%
    count(tracking_group) %>% filter(n >= 2) %>% pull(tracking_group)

  # resolve one display name per tracking group (Sapient by m/z -> putative -> m/z)
  lab <- sig %>% filter(tracking_group %in% incl) %>% group_by(tracking_group) %>%
    summarise(pn = { x <- putative_name[!is.na(putative_name) & putative_name != ""]; if (length(x)) x[1] else NA_character_ },
              mz2 = first(mz2), mz = first(mz), ion = first(ion_mode), .groups = "drop") %>%
    left_join(sap_lk, by = "mz2") %>%
    mutate(name = dplyr::coalesce(sap_name, pn),
           name = ifelse(is.na(name) | name %in% GENERIC_NAMES,
                         paste0("m/z ", sprintf("%.3f", mz), " (", ifelse(ion == "positive", "+", "-"), ")"),
                         name))

  # Per the Andrew x Trenton call (2026-08): show ONLY FDR-significant appearances.
  # One dot per compartment where the compound (tracking group) is significant.
  rep_cc <- sig %>% filter(tracking_group %in% incl) %>%
    group_by(tracking_group, compartment) %>% arrange(fdr) %>% slice(1) %>% ungroup()
  plot_df <- rep_cc %>%
    left_join(lab %>% select(tracking_group, row_id = name), by = "tracking_group") %>%
    group_by(row_id, compartment) %>% arrange(fdr) %>% slice(1) %>% ungroup()

  # Order compounds by the number of compartments reached (nsig), then mean effect.
  # The y factor is ascending, so ggplot places the 4-compartment ("headline")
  # transfers at the TOP, then 3-, then 2-compartment groups below.
  ord <- plot_df %>% group_by(row_id) %>%
    summarise(nsig = n(), meff = mean(effect_size), .groups = "drop") %>% arrange(nsig, meff)
  plot_df <- plot_df %>%
    mutate(row_id = factor(row_id, levels = ord$row_id),
           compartment = factor(compartment, levels = COMP_LEVELS))

  # Solid horizontal separators BETWEEN every component (each y tick, at k + 0.5);
  # the compartment-count group boundaries (4 -> 3 -> 2) are drawn heavier on top.
  # Gridlines are removed so these read cleanly.
  n_lv  <- nlevels(plot_df$row_id)
  row_y <- if (n_lv > 1) seq_len(n_lv - 1) + 0.5 else numeric(0)
  grp_y <- head(cumsum(rle(ord$nsig)$lengths), -1) + 0.5

  pd <- position_dodge(width = 0.7)
  p <- ggplot(plot_df, aes(effect_size, row_id, color = compartment, group = compartment)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    { if (length(row_y)) geom_hline(yintercept = row_y, color = "grey80", linewidth = 0.25) } +
    { if (length(grp_y)) geom_hline(yintercept = grp_y, color = "grey45", linewidth = 0.5) } +
    geom_errorbarh(aes(xmin = cil, xmax = ciu), position = pd, height = 0, linewidth = 0.4, alpha = 0.9) +
    geom_point(size = 2.1, alpha = 0.95, stroke = 0.7, position = pd) +
    scale_color_manual(values = COMP_COLS, name = "Compartment", drop = FALSE) +
    guides(colour = guide_legend(nrow = 2)) +
    labs(x = "Scaled average treatment effect (BEP), 95% CI", y = NULL) +
    theme_bw(base_size = 7) +
    theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 7), legend.text = element_text(size = 6),
          legend.title = element_text(size = 7), legend.key.size = unit(0.3, "cm"),
          panel.grid = element_blank(), legend.position = "bottom", legend.box = "vertical")

  # Right-column physical size (~92 mm) so 6-7 pt fonts print at size (no downscaling).
  ggsave(out_png, p, width = 106, height = 224, units = "mm", dpi = 300, device = ragg::agg_png)
  cat("wrote", out_png, "| compounds:", nlevels(plot_df$row_id), "\n")
  invisible(p)
}

if (sys.nframe() == 0) invisible(build_panelD())
