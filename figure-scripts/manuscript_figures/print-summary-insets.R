# print-summary-insets.R
# =============================================================================
# Two small inset panels for the print-summary figure (Trenton's layout):
#   (1) B-vitamins: density of milk concentrations, control vs intervention, for each
#       trial's most responsive B-vitamin at its strongest visit;
#   (2) Triglycerides: control -> BEP model-based mean concentration (uM) for the three
#       most significant triglycerides in each BEP trial at its strongest visit.
# Colours match Fig. 4 (control #0072B2, intervention #D55E00). Arm pooling matches the
# combined-arm analysis (src/2 analysis/2_adjusted_analysis_combined_arms.R):
#   ELICIT     Nico + Nico+Az.        vs Control + Az.
#   MISAME-III IFA/BEP + BEP/BEP      vs Control + BEP/IFA  (postnatal BEP)
#   Mumta-LW   BEP+ExBf + BEP+ExBf+AZT vs Control
# % labels = native-unit ATE / model-based control mean, the convention of the text
# (checked: ELICIT vitamin B3 at 1 mo gives the +183.7% quoted in the Results).
#
# Outputs: figures/print_summary_inset_bvitamins.{pdf,eps,png}
#          figures/print_summary_inset_triglycerides.{pdf,eps,png}
# Run from repo root: Rscript figure-scripts/manuscript_figures/print-summary-insets.R
# =============================================================================
suppressMessages({library(dplyr); library(ggplot2)})
root <- here::here()
source(file.path(root, "figure-scripts/0_figure-functions.R"))

ARM_COLS <- c(Control = "#0072B2", Intervention = "#D55E00")
inset_theme <- theme_imic(base_size = 8) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 8),
        plot.margin = margin(2, 4, 2, 2))

# ---- model-based % change (native ATE / model-based control mean) ----------------
nat <- readRDS(file.path(root, "results/adjusted_combined_arms_intervention_effects_unscaled_results_clean.RDS"))
pct_change <- function(st, vis, code) {
  x <- nat %>% filter(study == st, visit == vis, tolower(biomarker) == code)
  ctl <- x %>% filter(measure == "MN", contrast == "Control") %>% pull(est)
  ate <- x %>% filter(measure == "ATE") %>% pull(est)
  100 * ate[1] / ctl[1]
}

# ---- (1) B-vitamin densities ---------------------------------------------------------
d <- readRDS(file.path(root, "data/merged_analysis_datasets.RDS"))
INTERVENTION <- c("Nico", "Nico+Az.", "IFA/BEP", "BEP/BEP", "BEP+ExBf", "BEP+ExBf+AZT")
cells <- tibble::tribble(
  ~study,   ~visit, ~visit_lab, ~code, ~trial,       ~agent,
  "Elicit", "1",    "1 mo.",    "b3",  "ELICIT",     "nicotinamide",
  "Misame", "3",    "3-4 mo.",  "pa",  "MISAME-III", "BEP",
  "Vital",  "40",   "1.5 mo.",  "b2",  "Mumta-LW",   "BEP")
dens <- bind_rows(lapply(seq_len(nrow(cells)), function(i) {
  cc <- cells[i, ]
  d %>% filter(study == cc$study, as.character(visit) == cc$visit) %>%
    transmute(value = .data[[cc$code]],
              arm = factor(ifelse(arm %in% INTERVENTION, "Intervention", "Control"),
                           levels = names(ARM_COLS)),
              panel = i) %>%
    filter(!is.na(value), value > 0)
}))
cells$pct <- mapply(pct_change, cells$study, cells$visit_lab, cells$code)
print(cells %>% select(trial, code, visit_lab, pct))
# two-line plotmath strip: "ELICIT (nicotinamide)" / "Vitamin B3, 1 mo."
cells$strip <- sprintf('atop("%s (%s)", %s)', cells$trial, cells$agent,
                       plotmath_label(paste0(canonical_label(cells$code), ", ",
                                             sub("-", "–", cells$visit_lab))))
dens$strip <- factor(cells$strip[dens$panel], levels = cells$strip)
lab <- cells %>% mutate(strip = factor(strip, levels = cells$strip),
                        txt = sprintf("%+.0f%%", pct))
unit_x <- "Concentration, µg/L (log scale)"

p1 <- ggplot(dens, aes(x = value, fill = arm, colour = arm)) +
  geom_density(alpha = 0.45, linewidth = 0.4, adjust = 1.1) +
  geom_text(data = lab, aes(x = Inf, y = Inf, label = txt), inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.4, colour = ARM_COLS[["Intervention"]], fontface = "bold", size = 3) +
  facet_wrap(~ strip, nrow = 1, scales = "free", labeller = label_parsed) +
  scale_x_log10(breaks = function(lim) {            # 1-3-10 series: 2-4 ticks, no collisions
                  b <- 10^(floor(log10(lim[1])):ceiling(log10(lim[2]))); b <- sort(c(b, 3 * b))
                  b[b >= lim[1] & b <= lim[2]] },
                labels = scales::label_number(big.mark = ",", drop0trailing = TRUE)) +
  scale_fill_manual(values = ARM_COLS) + scale_colour_manual(values = ARM_COLS) +
  labs(x = unit_x, y = NULL) +
  inset_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_figure_3way(p1, "print_summary_inset_bvitamins", width = 5.2, height = 2.1,
                 dir = file.path(root, "figures"))

# ---- (2) triglyceride dumbbells ---------------------------------------------------------
tg <- read.csv(file.path(root, "results/tables/targeted_metabolites_native_units.csv")) %>%
  filter(grepl("^tg", tolower(biomarker)),
         (study == "Mumta-LW" & visit == "1.5 mo.") | (study == "MISAME-III" & visit == "3-4 mo."),
         qval_standardized < 0.05, ATE > 0) %>%
  group_by(study) %>% arrange(qval_standardized, .by_group = TRUE) %>% slice_head(n = 3) %>% ungroup()
tg_name <- function(b) {            # "Tg.18.3_36.2." -> "TG 54:5 (18:3)"
  m <- regmatches(b, regexec("^Tg\\.(\\d+)\\.(\\d+)_(\\d+)\\.(\\d+)\\.?$", b))
  vapply(m, function(v) if (length(v) == 5) {
    n <- as.integer(v[-1]); sprintf("TG %d:%d (%d:%d)", n[1] + n[3], n[2] + n[4], n[1], n[2])
  } else NA_character_, character(1))
}
tg <- tg %>% mutate(name = tg_name(biomarker),
                    trial = factor(sprintf("%s, %s", study, sub("-", "–", visit)),
                                   levels = c("Mumta-LW, 1.5 mo.", "MISAME-III, 3–4 mo.")),
                    name = factor(name, levels = rev(unique(name))),
                    txt = sprintf("%+.0f%%", pct_change))
print(tg %>% select(study, visit, biomarker, name, control_mean, intervention_mean, pct_change, qval_standardized))

p2 <- ggplot(tg, aes(y = name)) +
  geom_segment(aes(x = control_mean, xend = intervention_mean, yend = name), colour = "grey55",
               linewidth = 0.5, arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  geom_point(aes(x = control_mean, colour = "Control"), size = 2.2) +
  geom_point(aes(x = intervention_mean, colour = "Intervention"), size = 2.2) +
  geom_text(aes(x = intervention_mean, label = txt), hjust = -0.35, size = 2.6,
            colour = ARM_COLS[["Intervention"]], fontface = "bold") +
  facet_wrap(~ trial, ncol = 1, scales = "free") +
  scale_x_log10(labels = scales::label_number(big.mark = ",", drop0trailing = TRUE),
                expand = expansion(mult = c(0.08, 0.35))) +
  scale_colour_manual(values = ARM_COLS, labels = c(Control = "Control", Intervention = "BEP")) +
  labs(x = "Mean concentration, µM (log scale)", y = NULL) +
  inset_theme + theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
save_figure_3way(p2, "print_summary_inset_triglycerides", width = 3.0, height = 2.6,
                 dir = file.path(root, "figures"))
# ---- (3) treemap: which metabolite classes the interventions changed ------------------
# Tile area = number of FDR-significant targeted-metabolite effects (all trials and
# visits: 515, the Results' denominator; triglycerides = 73.4% of them). Colours are
# Fig. 5A's category colours. Drawn with plain rectangles (no treemap package needed):
# triglycerides fill the left block, the other classes stack in the right column.
eff <- read.csv(file.path(root, "results/tables/targeted_metabolites_native_units.csv")) %>%
  filter(qval_standardized < 0.05)
keep <- c("Triglycerides", "Diglycerides", "Acylcarnitines")
tm <- eff %>% mutate(class = ifelse(category %in% keep, category, "Other classes")) %>%
  count(class) %>% mutate(share = n / sum(n))
stopifnot(sum(tm$n) == 515)
TM_COLS <- c(Triglycerides = "#0072B2", Diglycerides = "#D55E00",
             Acylcarnitines = "#E69F00", `Other classes` = "#CC79A7")
p_tg <- tm$share[tm$class == "Triglycerides"]
right <- tm %>% filter(class != "Triglycerides") %>%
  mutate(class = factor(class, levels = c("Acylcarnitines", "Diglycerides", "Other classes"))) %>%
  arrange(class) %>% mutate(h = n / sum(n), ytop = 1 - cumsum(c(0, head(h, -1))), ybot = ytop - h)
rects <- bind_rows(
  tibble(class = "Triglycerides", xmin = 0, xmax = p_tg, ymin = 0, ymax = 1),
  right %>% transmute(class = as.character(class), xmin = p_tg, xmax = 1, ymin = ybot, ymax = ytop)) %>%
  left_join(tm, by = "class") %>%
  mutate(lab = sprintf("%s\n%d (%.0f%%)", sub(" classes", "\nclasses", class), n, 100 * share),
         x = (xmin + xmax) / 2, y = (ymin + ymax) / 2)
inside <- rects %>% filter(class != "Acylcarnitines")
ac <- rects %>% filter(class == "Acylcarnitines")
p3 <- ggplot(rects) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = class),
            colour = "white", linewidth = 0.8) +
  geom_text(data = inside, aes(x = x, y = y, label = lab,
                               colour = class == "Other classes"),
            size = ifelse(inside$class == "Triglycerides", 3.4, 2.5),
            fontface = ifelse(inside$class == "Triglycerides", "bold", "plain"), lineheight = 0.9) +
  # the acylcarnitine tile is a sliver (6 of 515): label it outside, with a leader
  geom_segment(data = ac, aes(x = 1.005, xend = 1.04, y = y, yend = y), colour = "grey40", linewidth = 0.3) +
  geom_text(data = ac, aes(x = 1.05, y = ymax, label = sprintf("Acylcarnitines\n%d (%.0f%%)", n, 100 * share)),
            hjust = 0, vjust = 1, size = 2.5, lineheight = 0.9) +   # hangs down from the tile top
  scale_fill_manual(values = TM_COLS, guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "grey15", `FALSE` = "white"), guide = "none") +
  scale_x_continuous(limits = c(0, 1.38), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(caption = sprintf("Tile area: FDR-significant effects on targeted milk metabolites (n = %d)", sum(tm$n))) +
  theme_void(base_size = 8) +
  theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey30"),
        plot.margin = margin(2, 2, 2, 2))
save_figure_3way(p3, "print_summary_inset_lipid_treemap", width = 3.2, height = 2.2,
                 dir = file.path(root, "figures"))
print(tm)
cat("wrote figures/print_summary_inset_{bvitamins,triglycerides,lipid_treemap}.{pdf,eps,png}\n")
