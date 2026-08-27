# =============================================================================
# 29-figure-overall-shift-agreement.R
#
# Figure C (overall directional agreement, beyond the joint-FDR features).
#   C1  among features BEP clearly moved in the MOTHER (p<0.05), the distribution
#       of the same features' INFANT effect, split by maternal direction -- the
#       maternal-up set sits positive, the maternal-down set negative.
#   C2  share of infant effects agreeing in direction as we tighten the maternal
#       side (all -> top 1000 -> top 200 -> maternal-FDR-sig). Two series:
#         same-run maternal<->infant blood (shared ID, no matching) -- rises;
#         cross-platform milk<->infant (25 ppm mass match) -- stays at chance
#         (negative control: agreement is recoverable only with exact linkage).
# Same-run pair = maternal & infant VAMS pn56 (combined-arms / any-postnatal-BEP).
# Within-dataset, direction-only. Out: figures/cross_compartment/
# =============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))
outdir <- paste0(root, "figures/cross_compartment/"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
up <- function(x) toupper(as.character(x))
TEAL <- "#1D9E75"; TEAL_D <- "#0F6E56"; GRAY <- "#888780"; GRAY_D <- "#5F5E5A"; INK <- "#2C2C2A"

milkC  <- as.data.table(readRDS(paste0(root,"results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")))
bloodC <- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))
# One row per milk feature at visit vv (keep the most significant duplicate):
# feature key, effect, p-value and FDR-significance flag.
gm <- function(vv){
  d <- milkC[measure=="ATE" & study=="Misame" & visit==vv]
  d[, feature := up(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]         # smallest-p row wins per feature
  d[, .(feature, est, pval, fdr=sigFDR)]
}
# Same shape for one blood dataset (dd) at one visit (vv).
gb <- function(dd, vv){
  d <- bloodC[measure=="ATE" & dataset==dd & visit==vv]
  d[, feature := up(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]         # smallest-p row wins per feature
  d[, .(feature, est, pval, fdr=sigFDR)]
}

# same-run pair: maternal pn56 <-> infant pn56 (shared vam_ id)
mm <- gb("VamsPostnatalMaternal","pn56"); ii <- gb("VamsPostnatalInfant","pn56")
same <- merge(mm[,.(feature,est_m=est,p_m=pval,fdr_m=fdr)], ii[,.(feature,est_i=est)], by="feature")
same <- same[is.finite(est_m)&is.finite(est_i)]

# cross-platform control: milk 1-2mo <-> infant pn12 (mass match, 25 ppm)
mk <- gm("1-2 mo."); inf12 <- gb("VamsPostnatalInfant","pn12")
sk <- best_match(mzrt_milk(misame_only=FALSE)[feature %in% mk$feature],
                 mzrt_vams()[feature %in% inf12$feature], use_rt=FALSE)
cross <- merge(sk, mk[,.(feature_A=feature,est_m=est,p_m=pval,fdr_m=fdr)], by="feature_A")
cross <- merge(cross, inf12[,.(feature_B=feature,est_i=est)], by="feature_B")
cross <- cross[is.finite(est_m)&is.finite(est_i)]

# ---- C2 stringency curve -----------------------------------------------------
# Directional agreement at four tightening cuts on the maternal side: all features,
# the top 1000 and top 200 by maternal p, and the maternal-FDR-significant set.
# (Named to avoid shadowing graphics::curve.)
stringency_curve <- function(d, lab){
  setorder(d, p_m)                                   # rank by maternal p so head() = strongest maternal effects
  agreement_pct <- function(s) if(nrow(s)>=5) round(100*mean(sign(s$est_m)==sign(s$est_i))) else NA_real_
  data.table(series=lab,
    stringency=factor(c("all","top 1k","top 200","mat. sig."), levels=c("all","top 1k","top 200","mat. sig.")),
    pct=c(agreement_pct(d), agreement_pct(head(d,1000)), agreement_pct(head(d,200)), agreement_pct(d[fdr_m==1])),
    n  =c(nrow(d), min(1000,nrow(d)), min(200,nrow(d)), sum(d$fdr_m==1)))
}
cc <- rbind(stringency_curve(same,"same run (exact match)"), stringency_curve(cross,"cross-platform (mass match)"))
cc[, series := factor(series, levels=c("same run (exact match)","cross-platform (mass match)"))]
cat("=== C2 stringency curve ===\n"); print(cc)

pC2 <- ggplot(cc, aes(stringency, pct, colour=series, group=series, linetype=series)) +
  geom_hline(yintercept=50, linetype="dotted", colour=GRAY) +
  annotate("text", x=Inf, y=52.5, label="chance", hjust=1.1, size=3, colour=GRAY_D) +
  geom_line(linewidth=1) + geom_point(size=2.6) +
  scale_colour_manual(values=c(TEAL, GRAY), name=NULL) +
  scale_linetype_manual(values=c("solid","42"), name=NULL) +
  scale_y_continuous(limits=c(45,100)) +
  labs(x="maternal-side stringency  →", y="infant agreement (% same direction)",
       title="C2   Agreement rises with effect size — only for exact linkage") +
  theme_minimal(base_size=12) +
  theme(legend.position=c(0.36,0.88), legend.background=element_rect(fill="white",colour=NA),
        legend.key.height=unit(11,"pt"), panel.grid.minor=element_blank(), plot.title=element_text(size=12.5))

# ---- C1 mirror density (same-run pair, maternal-moved features) --------------
dens <- copy(same[p_m<0.05]); dens[, mdir := factor(ifelse(est_m>0,"raised in mother","lowered in mother"),
                                                     levels=c("raised in mother","lowered in mother"))]
cat(sprintf("\nC1 density: maternal p<0.05 features = %d (up %d, down %d)\n",
            nrow(dens), sum(dens$est_m>0), sum(dens$est_m<0)))
pC1 <- ggplot(dens, aes(est_i, fill=mdir, colour=mdir)) +
  geom_vline(xintercept=0, linetype="dashed", colour=GRAY) +
  geom_density(alpha=0.35, linewidth=0.8) +
  scale_fill_manual(values=c(TEAL, GRAY), name=NULL) + scale_colour_manual(values=c(TEAL_D, GRAY_D), name=NULL) +
  coord_cartesian(xlim=c(-0.6, 0.85)) +
  labs(x="infant blood effect (SD)", y="density",
       title="C1   What BEP moved in the mother moves the same way in the infant") +
  theme_minimal(base_size=12) +
  theme(legend.position=c(0.74,0.88), legend.background=element_rect(fill="white",colour=NA),
        legend.key.height=unit(11,"pt"), panel.grid.minor=element_blank(), plot.title=element_text(size=12.5))

fig <- (pC1 | pC2) + plot_layout(widths=c(1,1)) +
  plot_annotation(title="Overall direction agrees across compartments where features are linked exactly",
                  subtitle="Same-run maternal↔infant blood, 5–6 mo. (combined-arms). Cross-platform pair shown as a negative control.",
                  theme=theme(plot.title=element_text(size=15), plot.subtitle=element_text(size=11, colour=GRAY_D)))
ggsave(paste0(outdir,"fig_overall_shift_agreement.png"), fig, width=12, height=5.6, dpi=200, bg="white")
cat("\nSaved figures/cross_compartment/fig_overall_shift_agreement.png\n")
