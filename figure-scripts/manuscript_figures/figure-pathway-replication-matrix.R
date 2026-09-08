# figure-pathway-replication-matrix.R
# =============================================================================
# Cross-study x cross-domain pathway-DIRECTION replication matrix.
#
# One tile per (pathway x study x analysis domain) that is FDR-significant
# (Q<0.05) in that cell, coloured by direction of the BEP effect in milk
# (blue = down, red = up) and shaded by -log10(FDR). Rows grouped into
# mechanistic THEMES; columns = analysis domain, x within-panel = study.
# The point is to make cross-study + cross-platform REPLICATION visible at a
# glance (a theme lit up in >1 column and >1 study is well supported).
#
# Metabolomics domains are FDR-significant enrichment cells collapsed to the
# most-significant timepoint/contrast per (pathway x study x domain).
# Proteome GO is PROVISIONAL (pre Trenton's protein-only rerun) and flagged.
#
# Out: figures/figure_pathway_replication_matrix.png
#      results/pathway_replication_matrix.csv
# =============================================================================
suppressMessages({ library(data.table); library(ggplot2) })
R <- paste0(here::here(), "/")
norm <- function(x) gsub("[^a-z0-9]+", "", tolower(x))   # lowercase FIRST (class is [a-z])

# --- pull FDR-significant enrichment cells from each metabolomics domain -------
grab <- function(path, pcol, scol, dcol, dom, fdrcol = "fdr_native") {
  d <- fread(path)
  if (!fdrcol %in% names(d) && "FDR" %in% names(d)) fdrcol <- "FDR"
  d <- d[get(fdrcol) < 0.05]
  if (!nrow(d)) return(NULL)
  data.table(domain = dom, study = as.character(d[[scol]]),
             direction = tolower(as.character(d[[dcol]])),
             pathway = as.character(d[[pcol]]), fdr = as.numeric(d[[fdrcol]]))
}
m <- rbindlist(list(
  grab(paste0(R,"results/metaboanalyst/primary_combined/primary_combined_supplementary_table.csv"),  "pathway","study","direction","Targeted\nprimary\n(B-vit.)"),
  grab(paste0(R,"results/metaboanalyst/tertiary_combined/tertiary_combined_supplementary_table.csv"), "pathway","study","direction","Targeted lipidome\n(tertiary)"),
  grab(paste0(R,"results/metaboanalyst/untargeted_msea/untargeted_msea_combined.csv"),                "pathway","study","direction","Untargeted\nMSEA"),
  grab(paste0(R,"results/metaboanalyst/mummichog_s5/milk_mummichog_tableS5.csv"),                     "Pathway","Study","Regulation","Untargeted\nMummichog")
), fill = TRUE)

# provisional proteome GO, keep a curated, representative subset of the S6 terms
# so the (large) proteome block does not visually swamp the metabolite themes.
pg <- fread(paste0(R,"results/metaboanalyst/proteomics_go/proteomics_go_tableS6.csv"))
pg_keep <- c("telomere maintenance", "RNA-templated DNA biosynthetic process",
             "DNA biosynthetic process", "translation", "cytoplasmic translation",
             "nucleotide biosynthetic process", "ribose phosphate biosynthetic process",
             "ribonucleoprotein complex biogenesis", "protein folding", "response to starvation",
             "negative regulation of peptidase activity", "vesicle-mediated transport")
pg <- pg[fdr < 0.05 & Description %in% pg_keep,
         .(domain = "Proteome GO\n(provisional)", study,
           direction = tolower(regulation), pathway = Description, fdr)]
m <- rbindlist(list(m, pg), fill = TRUE)

# --- normalise study labels + direction ---------------------------------------
m[, study := fifelse(grepl("Elicit|ELICIT", study), "ELICIT",
             fifelse(grepl("Misame|MISAME", study), "MISAME",
             fifelse(grepl("Vital|Mumta",  study), "Mumta", study)))]
m[, direction := fifelse(direction %in% c("down","downregulated"), "Down (BEP lowers)",
                 fifelse(direction %in% c("up","upregulated"),    "Up (BEP raises)", direction))]
m[, key := norm(pathway)]

# --- map pathways -> mechanistic themes + canonical labels --------------------
theme_of <- function(k) {
  first <- function(rules) { for (r in names(rules)) if (grepl(rules[[r]], k)) return(r); NA_character_ }
  first(c(
    "One-carbon / methyl donors"        = "methionine|glycineandserine|^glycine|betaine|homocysteine|folate|cysteine|serine",
    "Polyamine biosynthesis"            = "spermidine|spermine|polyamine",
    "NAD / nucleotide / PPP"            = "nicotin|purine|pyrimidine|pentosephosphate|nucleotide|ribosephosph|tryptophan",
    "Carnitine & fatty-acid oxidation"  = "carnitine|betaoxidation|branchedchainfatty|valineleucineisoleucine|fattyacidactivation|phytanic|saturatedfattyacids",
    "De novo lipogenesis & PUFA"        = "denovofattyacid|linoleate|linolenic|fattyacidelongation|fattyacidmetabolism",
    "Eicosanoid / inflammatory lipid"   = "arachidonic|leukotriene|prostaglandin",
    "Amino-acid & nitrogen handling"    = "arginineandproline|glutamate|glutamine|glutathione|ureacycle|ammonia|aspartate|histidine|tyrosine|betaalanine|lysine|threonine|selen",
    "Energy: TCA / glycolysis / ketone" = "citricacid|pyruvate|warburg|glycolysis|gluconeogenesis|ketonebody|butyrate|mitochondrialelectron|propanoate|heparansulfate",
    "Vitamin A / retinol"               = "retinol|vitamina",
    "Proteome: genome/telomere/transln" = "telomere|dnabiosynth|dnametabolic|translation|chromatin|chromosome|ribonucleoprotein|proteinfolding|proteinstability|proteinmaturation|nucleotidebiosynth|ribosephosphate|starvation|spermegg|cellrecognition|chaperone",
    "Proteome: immune / anti-proteolysis" = "immune|complement|endopeptidase|peptidase|hydrolase|proteolysis|humoral|vesicle"
  ))
}
m[, theme := vapply(key, theme_of, character(1))]
m <- m[!is.na(theme)]

# canonical short pathway label (merge a few synonym spellings)
m[, label := pathway]
m[grepl("denovofattyacid", key), label := "De novo fatty acid biosynthesis"]
m[grepl("nicotin", key) & domain %like% "primary", label := "Nicotinate & nicotinamide"]

# collapse to most-significant cell per (theme,label,study,domain)
setorder(m, fdr)
cell <- m[, .SD[1], by = .(theme, label, study, domain)]
cell[, neglog10fdr := pmin(-log10(fdr), 6)]

fwrite(cell[order(theme,label,domain,study),
            .(theme,label,study,domain=gsub("\n"," ",domain),direction,fdr,neglog10fdr)],
       paste0(R,"results/pathway_replication_matrix.csv"))

# --- order factors for the plot -----------------------------------------------
theme_order <- c("NAD / nucleotide / PPP","One-carbon / methyl donors","Polyamine biosynthesis",
                 "Amino-acid & nitrogen handling","Carnitine & fatty-acid oxidation",
                 "De novo lipogenesis & PUFA","Eicosanoid / inflammatory lipid",
                 "Energy: TCA / glycolysis / ketone","Vitamin A / retinol",
                 "Proteome: genome/telomere/transln","Proteome: immune / anti-proteolysis")
dom_order <- c("Targeted\nprimary\n(B-vit.)","Targeted lipidome\n(tertiary)",
               "Untargeted\nMSEA","Untargeted\nMummichog","Proteome GO\n(provisional)")
cell[, theme  := factor(theme, levels = rev(theme_order))]
cell[, domain := factor(domain, levels = dom_order)]
cell[, study  := factor(study, levels = c("ELICIT","MISAME","Mumta"))]
# row label unique within theme, ordered
lab_levels <- cell[order(theme, label), unique(label)]
cell[, label := factor(label, levels = rev(lab_levels))]

p <- ggplot(cell, aes(x = study, y = label)) +
  geom_tile(aes(fill = direction, alpha = neglog10fdr), color = "grey30", linewidth = 0.2) +
  facet_grid(theme ~ domain, scales = "free", space = "free", switch = "y") +
  scale_fill_manual(values = c("Down (BEP lowers)" = "#2C7FB8", "Up (BEP raises)" = "#D7301F"),
                    name = "Direction of BEP effect in milk") +
  scale_alpha_continuous(range = c(0.35, 1), name = expression(-log[10]*"(Q)"),
                         breaks = c(2,4,6), labels = c("2","4","≥6")) +
  labs(x = NULL, y = NULL,
       title = "BEP milk-metabolome pathway effects: cross-study × cross-platform replication",
       subtitle = "Tile = pathway FDR-significant (Q<0.05) in that study × domain; colour = direction, shade = significance") +
  theme_bw(base_size = 8) +
  theme(
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6.5, face = "bold"),
    strip.text.x = element_text(size = 6.8, face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    strip.placement = "outside",
    panel.spacing.y = unit(1.5, "pt"), panel.spacing.x = unit(2, "pt"),
    axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 6.6),
    panel.grid = element_blank(),
    plot.title = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 7),
    legend.position = "bottom", legend.box = "horizontal")

ggsave(paste0(R,"figures/figure_pathway_replication_matrix.png"), p,
       width = 9.5, height = 9, units = "in", dpi = 300, bg = "white")
cat("wrote figures/figure_pathway_replication_matrix.png |", nrow(cell), "significant cells across",
    uniqueN(cell$label), "pathways\n")
print(cell[, .N, by = .(domain)][order(domain)])
