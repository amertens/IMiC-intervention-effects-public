# =============================================================================
# appendix-blood-class-enrichment.R
# Direction-split chemical-class enrichment across the three blood compartments,
# as a dot plot to sit alongside the milk Fig 6. compartment (x) x class (y),
# faceted by direction (up/down); dot size = -log10(p), colour = direction,
# black ring = FDR < 0.05. Input: results/blood_chemical_class_enrichment_directional.csv
# (from 52-blood-class-enrichment-direction-sensitivity.R).
# =============================================================================
suppressMessages({library(data.table); library(ggplot2)})
root <- paste0(here::here(), "/")

d <- fread(paste0(root, "results/blood_chemical_class_enrichment_directional.csv"))
d[, neglogp := -log10(p)]
d[, fdr_sig := fdr < 0.05]
d[, comp := factor(compartment,
     levels = c("MaternalPlasma","VamsPostnatalMaternal","VamsPostnatalInfant","VamsPrenatal"),
     labels = c("Maternal\nplasma","Maternal\nVAMS","Infant\nVAMS","Prenatal VAMS\n(neg. control)"))]
d <- d[!is.na(comp)]
d[, direction := factor(direction, levels = c("up","down"), labels = c("up in BEP","down in BEP"))]

# keep classes with at least a nominal signal somewhere, order y by strongest
keep <- d[p < 0.10, unique(class)]
d <- d[class %in% keep]
ord <- d[, .(m = max(neglogp)), by = class][order(m)]$class
d[, class := factor(class, levels = ord)]

p <- ggplot(d, aes(comp, class)) +
  geom_point(aes(size = neglogp, fill = direction, colour = fdr_sig), shape = 21, stroke = 0.9) +
  facet_wrap(~direction) +
  scale_fill_manual(values = c("up in BEP" = "#B23A48", "down in BEP" = "#2C6E9E"), guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75"), name = "FDR < 0.05") +
  scale_size_continuous(range = c(1.5, 9), name = expression(-log[10](p))) +
  labs(x = NULL, y = NULL,
       title = "BEP-responsive chemical-class enrichment in blood (direction-split)",
       subtitle = paste("Fisher over-representation vs all classed features. Black ring = FDR<0.05.",
                        "Plasma classes are m/z-borrowed (isobaric); VAMS classes are direct.", sep = "\n")) +
  theme_bw(base_size = 10) +
  theme(plot.subtitle = element_text(colour = "#666", size = 8),
        panel.grid.minor = element_blank(), legend.position = "right")

# Write the published supplementary-figure name (Fig S8) directly, plus a working copy in results/.
ggsave(paste0(root, "figures/appendix_blood_class_enrichment.png"), p, width = 10, height = 6.2, dpi = 200, bg = "white")
ggsave(paste0(root, "results/blood_class_enrichment_figure.png"),   p, width = 10, height = 6.2, dpi = 200, bg = "white")
cat("saved figures/appendix_blood_class_enrichment.png (published Fig S8) + results/ working copy (", length(keep), "classes)\n")
