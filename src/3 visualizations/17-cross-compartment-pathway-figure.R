# =============================================================================
# 17-cross-compartment-pathway-figure.R
#
# Cross-compartment Mummichog pathway comparison: MILK vs MATERNAL BLOOD vs INFANT
# BLOOD. Reads every mcg_pathwayanalysis_*.xlsx produced by the milk + blood
# mummichog runs (scripts 15/16), collapses sources to three compartment groups,
# and highlights pathways enriched across compartments -- the supplement figure
# for the milk -> maternal-blood -> infant-blood transfer story.
#
# Groups: Milk | Maternal blood (plasma + maternal VAMS) | Infant blood (infant VAMS).
# Outputs (tag-aware via BLOOD_ADJUST):
#   figures/cross_compartment/pathway_milk_vs_blood<_adjusted>.png
#   results/cross_compartment_pathways<_adjusted>.csv
# =============================================================================

suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(readxl)})
root <- paste0(here::here(), "/")
if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- FALSE
.osuf <- if (BLOOD_ADJUST) "_adjusted" else ""
outfig <- paste0(root, "figures/cross_compartment"); dir.create(outfig, showWarnings = FALSE, recursive = TRUE)

# milk pathways always come from the (unadjusted-dir) milk run; blood from the
# matching adjusted/unadjusted mummichog output dir.
dirs <- unique(c(paste0(root, "results/mummichog_output"),
                 paste0(root, "results/mummichog_output", .osuf)))
fs <- unlist(lapply(dirs, function(d) list.files(d, pattern = "mcg_pathwayanalysis_.*\\.xlsx$",
                                                 recursive = TRUE, full.names = TRUE)))
stopifnot(length(fs) > 0)

all <- bind_rows(lapply(fs, function(f) {
  d <- suppressMessages(read_excel(f))
  d$source <- sub("^[0-9.]+\\.", "", basename(dirname(dirname(f))))
  d
})) %>%
  rename(p = `p-value`) %>%
  mutate(group = case_when(
    grepl("^Milk", source)                               ~ "Milk",
    grepl("MaternalPlasma|VamsPostnatalMaternal", source)~ "Maternal blood",
    grepl("VamsPostnatalInfant", source)                 ~ "Infant blood",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(group))

# best (smallest) p per pathway per compartment group
g <- all %>% group_by(pathway, group) %>%
  summarise(p = min(p, na.rm = TRUE),
            enr = max(overlap_size / pathway_size, na.rm = TRUE), .groups = "drop") %>%
  mutate(sig = p < 0.05, neglog10p = pmin(-log10(p), 6))

# wide table + cross-compartment count
wide <- g %>% select(pathway, group, p) %>% pivot_wider(names_from = group, values_from = p)
ncomp <- g %>% group_by(pathway) %>% summarise(n_sig_groups = sum(sig),
                                               milk_sig = any(sig[group == "Milk"]), .groups = "drop")
tab <- left_join(wide, ncomp, by = "pathway") %>% arrange(desc(n_sig_groups))
write.csv(tab, paste0(root, "results/cross_compartment_pathways", .osuf, ".csv"), row.names = FALSE)

# figure: pathways enriched in >=2 compartments (the shared, transfer-relevant set)
shared <- ncomp %>% filter(n_sig_groups >= 2) %>% pull(pathway)
gd <- g %>% filter(pathway %in% shared) %>%
  mutate(group = factor(group, levels = c("Milk", "Maternal blood", "Infant blood")),
         pathway = factor(pathway, levels = rev(ncomp %>% filter(n_sig_groups >= 2) %>%
                                                  arrange(n_sig_groups, milk_sig) %>% pull(pathway))))
p <- ggplot(gd, aes(group, pathway, fill = neglog10p)) +
  geom_tile(colour = "grey90") +
  geom_text(data = filter(gd, sig), aes(label = "*"), size = 5, vjust = 0.78) +
  scale_fill_gradient(low = "white", high = "#E69F00", name = "-log10 p\n(capped 6)") +
  labs(x = NULL, y = NULL,
       title = paste0("Cross-compartment pathway enrichment (BEP)", if (BLOOD_ADJUST) " - adjusted" else ""),
       subtitle = "* = enriched (p<0.05); shown: pathways enriched in >=2 compartments") +
  theme_minimal(base_size = 10) + theme(axis.text.y = element_text(size = 8),
                                        panel.grid = element_blank())
ggsave(p, file = paste0(outfig, "/pathway_milk_vs_blood", .osuf, ".png"),
       width = 7.5, height = 0.28 * length(shared) + 2, limitsize = FALSE)

cat("Cross-compartment pathways: ", length(shared), " enriched in >=2 compartments",
    " (of ", dplyr::n_distinct(g$pathway), " total).\n", sep = "")
cat("Wrote figures/cross_compartment/pathway_milk_vs_blood", .osuf, ".png and ",
    "results/cross_compartment_pathways", .osuf, ".csv\n", sep = "")
