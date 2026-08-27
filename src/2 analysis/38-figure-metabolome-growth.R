# =============================================================================
# 38-figure-metabolome-growth.R   Figure H: does the BEP-responsive infant
# carnitine (vam_1005523, pn12) relate to the infant's later growth? Exploratory,
# Kim-independent. Links infant carnitine (idBiospe) -> subjid -> growth z-scores.
# =============================================================================
suppressMessages({library(data.table); library(ggplot2); library(haven)})
root <- paste0(here::here(), "/"); add <- paste0(root,"data/additional datasets/"); out <- paste0(root,"figures/cross_compartment/")
TEAL<-"#1D9E75"; TEALD<-"#0F6E56"; GRAY<-"#888780"

# infant carnitine pn12 per idBiospe
vm <- fread(paste0(root,"results/_carn_vam.csv")); setnames(vm,1,"feature")
vl <- melt(vm, id.vars="feature", variable.name="s", value.name="carnitine"); vl[, s:=as.character(s)]
vp <- tstrsplit(vl$s,"_",fixed=TRUE); vl[, `:=`(v=vp[[1]], idBiospe=suppressWarnings(as.numeric(vp[[3]])), dy=fifelse(vp[[4]]=="e","infant","mother"))]
inf <- vl[v=="pn12" & dy=="infant" & is.finite(idBiospe), .(idBiospe, carnitine)][!duplicated(idBiospe)]

# idBiospe -> subjid bridge (from milk merged data; the trailing number of bmid,
# after the last "_" or "-", is the idBiospe that keys the carnitine samples).
mm <- as.data.table(readRDS(paste0(root,"data/merged_analysis_datasets.RDS")))
br <- unique(mm[study=="Misame", .(idBiospe=as.numeric(sub(".*[_-]","",bmid)), subjid=as.integer(subjid))])[!is.na(idBiospe)&!is.na(subjid)]

# growth (MISAME-3), child z-scores
g <- fread(paste0(root,"results/imic_growth_outcomes_dataset.csv"))
g <- g[studyid=="MISAME-3", .(subjid=as.integer(subjid), whz_6mo, waz_6mo, haz_6mo,
                              wt_vel=wt_vel_z_birth_6mo, len_vel=len_vel_z_birth_6mo)]

d <- merge(merge(inf, br, by="idBiospe"), g, by="subjid")
outcomes <- c(whz_6mo="weight-for-length z (6 mo)", waz_6mo="weight-for-age z (6 mo)",
              haz_6mo="length-for-age z (6 mo)", wt_vel="weight velocity z (birth-6 mo)",
              len_vel="length velocity z (birth-6 mo)")
long <- melt(d, id.vars=c("subjid","carnitine"), measure.vars=names(outcomes),
             variable.name="outcome", value.name="z")[is.finite(z)&is.finite(carnitine)]
long[, lab := outcomes[as.character(outcome)]]
stats <- long[, { ct<-cor.test(carnitine,z); .(n=.N, r=round(unname(ct$estimate),2), p=signif(ct$p.value,2)) }, by=lab]
cat("=== infant carnitine (pn12) vs later growth ===\n"); print(stats)
labstat <- merge(long[, .(x=min(carnitine), y=max(z)), by=lab], stats, by="lab")
labstat[, txt := sprintf("r=%.2f, p=%.2g (n=%d)", r, p, n)]

pH <- ggplot(long, aes(carnitine, z)) +
  geom_point(colour=GRAY, alpha=0.5, size=1.2) +
  geom_smooth(method="lm", se=TRUE, colour=TEALD, fill=TEAL, alpha=0.15, linewidth=0.8) +
  geom_text(data=labstat, aes(x=x, y=y, label=txt), hjust=0, vjust=1, size=3, colour="#2C2C2A") +
  facet_wrap(~lab, scales="free", nrow=2) +
  labs(x="infant blood carnitine at 1-2 mo (standardized)", y="growth z-score",
       title="H   Infant carnitine at 1-2 mo vs later growth (exploratory)",
       subtitle="MISAME-3 infants. Association of the BEP-responsive carnitine with anthropometry; correlations shown per panel.") +
  theme_minimal(base_size=10) + theme(plot.title=element_text(size=12),
       plot.subtitle=element_text(size=9.5, colour=GRAY), strip.text=element_text(size=9))
ggsave(paste0(out,"fig_H_metabolome_growth.png"), pH, width=9, height=6, dpi=200, bg="white")
fwrite(stats, paste0(root,"results/infant_carnitine_growth_assoc.csv"))
cat(sprintf("\nSaved fig_H_metabolome_growth.png (n dyads with growth = %d)\n", nrow(d)))
