# =============================================================================
# 26-figure-transfer-story.R
#
# Main-manuscript candidate figure (Framing A): the milk-mediated transfer story.
#   Panel A  schematic of the randomized BEP -> maternal blood -> milk -> infant path
#   Panel B  m/z 286.202 medium-chain acylcarnitine effect (combined-arms /
#            any-postnatal-BEP, 95% CI). The compound is deliberately NOT named:
#            the assignment is mass-based with no MS/MS confirmation, and the
#            manuscript reports only the directional cross-compartment increase.
#            across milk, maternal plasma, maternal blood-microsample, infant blood
#   Panel C  count of FDR-significant INCREASED infant features over the breastfeeding
#            visits, by arm -- grows only in the two post-birth-BEP arms
#
# Effects are within-dataset, standardized; panel B compares DIRECTION + significance,
# not absolute magnitude (Kim constraint). Out: figures/cross_compartment/
# =============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
root <- paste0(here::here(), "/")
outdir <- paste0(root, "figures/cross_compartment/"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
low <- function(x) tolower(as.character(x))
# DECK_VARIANT=1 -> drop the descriptive panel titles (keep A/B/C letters) for slides
DECK <- nzchar(Sys.getenv("DECK_VARIANT"))

milkC  <- as.data.table(readRDS(paste0(root,"results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")))
bloodC <- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))
blood  <- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_intervention_effects_results_clean.RDS")))

TEAL <- "#1D9E75"; TEAL_D <- "#0F6E56"; TEAL_L <- "#9FE1CB"; GRAY <- "#888780"; GRAY_L <- "#D3D1C7"; INK <- "#2C2C2A"

# ---- Panel B: carnitine across compartments (combined-arms primary) ----------
# One row per compartment for the SAME carnitine molecule. Each platform catalogues
# it under its own feature id (rlc_pos_mtb_* in milk/plasma, vam_* in the VAMS runs),
# so the ids differ by row; we pull each from the matching dataset/visit.
pb <- rbindlist(list(
  data.table(compartment="Breast milk\n(1–2 mo.)",            milkC[measure=="ATE" & study=="Misame" & low(biomarker)=="rlc_pos_mtb_3033006" & visit=="1-2 mo.", .(est,cil,ciu)]),
  data.table(compartment="Maternal plasma\n(1–2 mo.)",        bloodC[measure=="ATE" & dataset=="MaternalPlasma" & low(biomarker)=="rlc_pos_mtb_2679130" & visit=="pn12", .(est,cil,ciu)]),
  data.table(compartment="Maternal blood\nmicrosample (5–6 mo.)", bloodC[measure=="ATE" & dataset=="VamsPostnatalMaternal" & low(biomarker)=="vam_1005523" & visit=="pn56", .(est,cil,ciu)]),
  data.table(compartment="Infant blood\nmicrosample (1–2 mo.)",   bloodC[measure=="ATE" & dataset=="VamsPostnatalInfant" & low(biomarker)=="vam_1005523" & visit=="pn12", .(est,cil,ciu)])
))
pb[, compartment := factor(compartment, levels=rev(compartment))]

pB <- ggplot(pb, aes(est, compartment)) +
  geom_vline(xintercept=0, linetype="dashed", colour=GRAY) +
  geom_errorbarh(aes(xmin=cil, xmax=ciu), height=0.18, colour=TEAL_D, linewidth=0.7) +
  geom_point(size=3.2, colour=TEAL) +
  geom_text(aes(x=ciu, label=sprintf("+%.2f", est)), hjust=-0.25, size=3.5, colour=INK) +
  scale_x_continuous(limits=c(0, 1.85), breaks=seq(0,1.5,0.5)) +
  labs(x="BEP effect, standardized (SD units, 95% CI)", y=NULL,
       title=if(DECK) "B" else "B   Increased in every compartment (m/z 286.202)") +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        plot.title=element_text(size=12.5, face="plain"), axis.text.y=element_text(colour=INK, lineheight=0.95))

# ---- Panel C: infant FDR-sig INCREASED features over time, by arm ------------
inf <- blood[measure=="ATE" & dataset=="VamsPostnatalInfant" & sigFDR==1 & visit %in% c("pn12","pn34","pn56") & est>0]
cnt <- inf[, .(n_up=.N), by=.(visit, contrast)]
# Complete the visit x arm grid so combinations with no significant features plot
# as an explicit 0 rather than a missing point (keeps each arm's line continuous).
grid <- CJ(visit=c("pn12","pn34","pn56"), contrast=c("BEP/BEP","IFA/BEP","BEP/IFA"))
cnt <- merge(grid, cnt, by=c("visit","contrast"), all.x=TRUE); cnt[is.na(n_up), n_up:=0L]
arm_map <- c("BEP/BEP"="BEP both periods", "IFA/BEP"="BEP after birth only", "BEP/IFA"="BEP in pregnancy only")
cnt[, arm := factor(arm_map[contrast], levels=arm_map)]
cnt[, vis := factor(c(pn12="1–2", pn34="3–4", pn56="5–6 mo.")[visit], levels=c("1–2","3–4","5–6 mo."))]
arm_cols <- c("BEP both periods"=TEAL_D, "BEP after birth only"=TEAL, "BEP in pregnancy only"=GRAY)

pC <- ggplot(cnt, aes(vis, n_up, colour=arm, group=arm)) +
  geom_line(aes(linetype=arm), linewidth=0.9) +
  geom_point(size=2.8) +
  scale_colour_manual(values=arm_cols, name=NULL) +
  scale_linetype_manual(values=c("solid","solid","31"), name=NULL) +
  labs(x="infant visit", y="increased features (FDR<0.05)",
       title=if(DECK) "C" else "C   Infant signal grows over breastfeeding (post-birth arms)") +
  theme_minimal(base_size=12) +
  theme(legend.position=c(0.30,0.82), legend.background=element_rect(fill="white", colour=NA),
        legend.key.height=unit(11,"pt"), panel.grid.minor=element_blank(),
        plot.title=element_text(size=12.5))

# ---- Panel A: path schematic -------------------------------------------------
bx <- data.table(
  lab=c("Maternal\nBEP","Maternal\nblood","Breast\nmilk","Infant\nblood"),
  cx=c(1.4,4.2,7.0,9.8), fill=c(GRAY_L, TEAL_L, TEAL_L, TEAL_L))
bw <- 1.05; bh <- 0.62
pA <- ggplot(bx) +
  geom_rect(aes(xmin=cx-bw, xmax=cx+bw, ymin=1-bh, ymax=1+bh, fill=I(fill)), colour=GRAY) +
  geom_text(aes(cx, 1, label=lab), colour=INK, size=4, lineheight=0.9) +
  geom_segment(data=data.table(x=c(2.55,5.35,8.15)),
               aes(x=x, xend=x+0.85, y=1, yend=1), colour=TEAL_D, linewidth=0.9,
               arrow=arrow(length=unit(7,"pt"), type="closed")) +
  (if (!DECK) annotate("text", x=5.6, y=1.95, label="randomized supplement · transfer carried by the two post-birth-BEP arms",
           size=3.6, colour=GRAY, fontface="italic")) +
  scale_x_continuous(limits=c(0,11.2)) + scale_y_continuous(limits=c(0.2,2.1)) +
  labs(title=if(DECK) "A" else "A   One molecule, one path") +
  theme_void(base_size=12) + theme(plot.title=element_text(size=12.5, hjust=0, margin=margin(b=2)))

fig <- pA / (pB | pC) + plot_layout(heights=c(0.75, 2)) +
  plot_annotation(title="BEP transfers along maternal blood → milk → infant blood",
                  theme=theme(plot.title=element_text(size=15)))

out_fn <- if(DECK) "fig_transfer_story_deck.png" else "fig_transfer_story.png"
ggsave(paste0(outdir,out_fn), fig, width=11, height=7.2, dpi=200, bg="white")
cat("Saved figures/cross_compartment/", out_fn, "\n")
print(pb); print(cnt[order(arm,vis)])
