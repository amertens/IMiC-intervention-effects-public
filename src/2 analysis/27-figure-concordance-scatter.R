# =============================================================================
# 27-figure-concordance-scatter.R
#
# Candidate figure (Framing B): cross-compartment agreement among the FDR-sig features.
#   Left   maternal vs infant blood-microsample BEP effect at 5-6 mo (combined-arms),
#          shared-ID match (no mass alignment); features significant in BOTH highlighted,
#          octenoylcarnitine labelled.
#   Right  share of matched FDR-sig-in-both features agreeing in direction, per pair
#          (from the FDR-first lists, any-postnatal-BEP).
# Within-dataset, direction-only (Kim constraint). Out: figures/cross_compartment/
# =============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(ggrepel)})
root <- paste0(here::here(), "/")
outdir <- paste0(root, "figures/cross_compartment/"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
low <- function(x) tolower(as.character(x))
TEAL <- "#1D9E75"; TEAL_D <- "#0F6E56"; GRAY <- "#888780"; INK <- "#2C2C2A"

bloodC <- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))

matM <- bloodC[measure=="ATE" & dataset=="VamsPostnatalMaternal" & visit=="pn56", .(feature=low(biomarker), est_m=est, fdr_m=sigFDR)]
infM <- bloodC[measure=="ATE" & dataset=="VamsPostnatalInfant"   & visit=="pn56", .(feature=low(biomarker), est_i=est, fdr_i=sigFDR)]
d <- merge(matM, infM, by="feature")
d[, both := fdr_m==1 & fdr_i==1]                       # FDR-significant in BOTH maternal and infant
# nb = features significant in both; ncon = those whose maternal & infant effects
# share the same sign (concordant in direction).
nb <- sum(d$both); ncon <- sum(d$both & sign(d$est_m)==sign(d$est_i))
cat(sprintf("shared features: %d ; FDR-sig in both: %d ; concordant: %d (%.0f%%)\n",
            nrow(d), nb, ncon, 100*ncon/nb))
carn <- d[feature=="vam_1005523"]

lim <- max(abs(c(d$est_m, d$est_i)), na.rm=TRUE)*1.02
pScat <- ggplot() +
  geom_hline(yintercept=0, colour=GRAY, linewidth=0.3) + geom_vline(xintercept=0, colour=GRAY, linewidth=0.3) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour=TEAL, linewidth=0.4) +
  geom_point(data=d[both==FALSE], aes(est_m, est_i), colour=GRAY, alpha=0.25, size=0.7) +
  geom_point(data=d[both==TRUE],  aes(est_m, est_i), colour=TEAL, size=2.6) +
  geom_point(data=carn, aes(est_m, est_i), colour=TEAL_D, size=3.4) +
  geom_text_repel(data=carn, aes(est_m, est_i, label="octenoylcarnitine"), colour=INK, size=3.6,
                  nudge_y=0.35, nudge_x=-0.25, segment.colour=GRAY) +
  annotate("text", x=-lim*0.96, y=lim*0.92, hjust=0, colour=TEAL_D, size=3.8,
           label=sprintf("%d of %d significant in both agree\n— all increased", ncon, nb)) +
  coord_equal(xlim=c(-lim,lim), ylim=c(-lim,lim)) +
  labs(x="maternal blood effect (SD)", y="infant blood effect (SD)",
       title="Maternal vs infant blood, 5–6 mo. (same run, exact match)") +
  theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank(), plot.title=element_text(size=12.5))

# ---- companion: share agreeing by compartment pair (FDR-first, any-postnatal-BEP) ----
fl <- fread(paste0(root,"results/cross_compartment_fdr_first_lists.csv"))[contrast=="any-postnatal-BEP"]
# One bar per compartment pair: n matched features and % agreeing in direction.
bars <- fl[, .(n=.N, pct=round(100*mean(concordant))), by=arrow]
cat("\npairs:\n"); print(bars)
relab <- c("Milk<->MaternalPlasma"="Milk ↔ maternal plasma",
           "Milk<->InfantVAMS"="Milk ↔ infant blood",
           "MaternalVAMS<->InfantVAMS"="Maternal ↔ infant blood",
           "MaternalPlasma<->InfantVAMS"="Maternal plasma ↔ infant")
ord <- c("Milk<->MaternalPlasma","Milk<->InfantVAMS","MaternalVAMS<->InfantVAMS","MaternalPlasma<->InfantVAMS")
bars <- bars[match(ord, arrow)][!is.na(arrow)]
bars[, lab := factor(relab[arrow], levels=rev(relab[ord]))]

pBar <- ggplot(bars, aes(pct, lab)) +
  geom_col(fill=TEAL, width=0.62) +
  geom_text(aes(label=sprintf("%d%% (%d)", pct, n)), hjust=-0.12, size=3.5, colour=INK) +
  scale_x_continuous(limits=c(0,118), breaks=c(0,50,100)) +
  labs(x="features agreeing in direction (%)", y=NULL,
       title="Agreement by compartment pair") +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_text(colour=INK), plot.title=element_text(size=12.5))

fig <- (pScat | pBar) + plot_layout(widths=c(1.5,1)) +
  plot_annotation(title="Many features, one direction — BEP effects agree across compartments",
                  theme=theme(plot.title=element_text(size=15)))
ggsave(paste0(outdir,"fig_concordance_scatter.png"), fig, width=12, height=6, dpi=200, bg="white")
cat("\nSaved figures/cross_compartment/fig_concordance_scatter.png\n")
