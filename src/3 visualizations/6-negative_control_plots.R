# =============================================================================
# src/3 visualizations/6-negative_control_plots.R
#
# Reads:  results/unadjusted_intervention_effects_results.RDS
#         results/unadjusted_intervention_effects_results_blinded.RDS
# Writes: results/negative_control_results.RDS
#
# Paths above were recovered from this script's syntax tree and are
# repo-relative; they resolve from the repo root via here::here().
#
# Header generated from the code itself; it makes no claim about method.
# See README.md for run order and results/ARTIFACT_MANIFEST.csv for the
# exhibit each script feeds.
# =============================================================================


rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

d <- readRDS(paste0(here::here(),"/results/unadjusted_intervention_effects_results.RDS"))

#Use Elicit antibiotics as negative control (they were given to infants after BM collection) and use 
nc_res1 <- d$res_primary$res$`Elicit-1`$res %>% filter(contrast=="Az.", measure=="ATE")
nc_res2 <- d$res_secondary$res$`Elicit-1`$res %>% filter(contrast=="Az.", measure=="ATE")
nc_res3 <- d$res_tertiary$res$`Elicit-1`$res %>% filter(contrast=="Az.", measure=="ATE")
nc_res <- bind_rows(nc_res1,nc_res2,nc_res3)

#mark sig treatment effect
nc_res <- nc_res %>% mutate(sig_int = 1*(cil<0 & ciu<0 | cil>0 & ciu>0))
dim(nc_res)
as.numeric(prop.table(table(nc_res$sig_int))[2]*100)

# b3 associated components and catabolites
nico_vars = c("b3","nam","nmn","nr","nad","nufa")

#Elicit Nico as positive control (Vitamin b supplementation that should have an effect specifically on the B-vitamins).
pc_res1 <- d$res_primary$res$`Elicit-1`$res %>% 
  filter(contrast=="Nico", measure=="ATE", biomarker %in% nico_vars) %>%
  mutate(Visit="1 month")
pc_res5 <- d$res_primary$res$`Elicit-1`$res %>% 
  filter(contrast=="Nico", measure=="ATE", biomarker %in% nico_vars) %>%
  mutate(Visit="5 month")
pc_res <- bind_rows(pc_res1,pc_res5)

elicit_pc_plot <- ggplot(pc_res, 
                             aes(x=biomarker, y=est)) + geom_point() +
  geom_linerange(aes(ymin=cil, ymax=ciu)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  coord_flip() +
  facet_wrap(~Visit, scale="free_y", ncol=3) +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6),
        legend.position="none")  + xlab("B vitamin") + ylab("Z-score difference") + theme_bw()


saveRDS(list(nc_res=nc_res, 
             pc_res=pc_res,
             elicit_pc_plot=elicit_pc_plot),
        file=paste0(here::here(),"/results/negative_control_results.RDS"))


#check that blinded significance is ~5%
d_blinded <- readRDS(paste0(here::here(),"/results/unadjusted_intervention_effects_results_blinded.RDS"))

blinded_res1 <- do.call(rbind.data.frame, lapply(d_blinded$res_primary$res, `[[`, "res")) %>% filter(measure=="ATE")
blinded_res2 <- do.call(rbind.data.frame, lapply(d_blinded$res_secondary$res, `[[`, "res")) %>% filter(measure=="ATE")
blinded_res3 <- do.call(rbind.data.frame, lapply(d_blinded$res_tertiary$res, `[[`, "res")) %>% filter(measure=="ATE")
blinded_res <- bind_rows(blinded_res1,blinded_res2,blinded_res3) %>% mutate(sig_int = 1*(cil<0 & ciu<0 | cil>0 & ciu>0))
as.numeric(prop.table(table(blinded_res$sig_int))[2]*100)


d_unblinded <- readRDS(paste0(here::here(),"/results/unadjusted_intervention_effects_results.RDS"))

unblinded_res1 <- do.call(rbind.data.frame, lapply(d_unblinded$res_primary$res, `[[`, "res")) %>% filter(measure=="ATE")
unblinded_res2 <- do.call(rbind.data.frame, lapply(d_unblinded$res_secondary$res, `[[`, "res")) %>% filter(measure=="ATE")
unblinded_res3 <- do.call(rbind.data.frame, lapply(d_unblinded$res_tertiary$res, `[[`, "res")) %>% filter(measure=="ATE")
unblinded_res <- bind_rows(unblinded_res1,unblinded_res2,unblinded_res3) %>% mutate(sig_int = 1*(cil<0 & ciu<0 | cil>0 & ciu>0))
as.numeric(prop.table(table(unblinded_res$sig_int))[2]*100)
