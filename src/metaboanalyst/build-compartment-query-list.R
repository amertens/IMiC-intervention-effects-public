# =============================================================================
# build-compartment-query-list.R
#
# Builds the compound query list for the cross-compartment KEGG pathway analysis
# (Table S11). Starting from the FDR-significant MISAME-III features in maternal
# blood, milk and infant blood, it keeps up-regulated features that are linked to
# a same-direction feature in another compartment (ct_compartment_pairs(): shared
# feature ID, shared putative name, or same ion mode within 25 ppm), groups them
# by their own putative name (ct_name_track()), and keeps names seen in at least
# two compartment x timepoint cells. The console summary gives the counts behind
# the stated result (20 up-regulated features across at least two compartments,
# 16 of them reliably detected in the BEP supplement, 2 in all three).
# Run this before src/metaboanalyst/run-compartment-pathway.R.
#
# Inputs : results/supplement_status_fdr_features.csv  (src/2 analysis/54-supplement-detection-status.R)
#          src/metaboanalyst/R/compartment-tracking.R
# Outputs: results/compartment_tracking/linked_upregulated.csv    (row-level data, one row per name x cell)
#          results/compartment_tracking/pathway_compound_list.csv (Table S11 query list)
#
# Run from repo root: Rscript src/metaboanalyst/build-compartment-query-list.R
# =============================================================================
suppressMessages({ library(dplyr); library(tidyr); library(readr) })
source("src/metaboanalyst/R/compartment-tracking.R")   # ct_compartment_pairs(), ct_name_track()

SIGFEAT  <- "results/supplement_status_fdr_features.csv"
ROWDATA  <- "results/compartment_tracking/linked_upregulated.csv"
CMPDLIST <- "results/compartment_tracking/pathway_compound_list.csv"

# compartment x timepoint order used in the original R Markdown analysis
CT_LEVELS <- c("Maternal plasma · incl","Maternal plasma · tri3","Maternal plasma · pn12",
               "Maternal VAMS · tri3","Maternal VAMS · pn56","Milk · 1421d","Milk · pn12",
               "Milk · pn34","Infant VAMS · acco","Infant VAMS · pn12","Infant VAMS · pn34",
               "Infant VAMS · pn56")
# Maternal plasma and maternal VAMS are two assays of one compartment, so the
# "across compartments" counts use maternal blood / milk / infant blood.
POOL3 <- c("Maternal plasma" = "Maternal blood", "Maternal VAMS" = "Maternal blood",
           "Milk" = "Milk", "Infant VAMS" = "Infant blood")
# A bare compound class carried by the upstream name map; not a compound identity,
# so it is left out of the pathway query.
JUNK_NAMES <- c("Acid")

is_true <- function(x) { if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes") }

sig <- read_csv(SIGFEAT, show_col_types = FALSE) %>% filter(!is.na(mz), !is.na(effect_size))

# features that have at least one same-direction partner in another compartment
pairs   <- ct_compartment_pairs(sig)
tracked <- pairs %>% filter(same_direction) %>%
  select(tracking_id_1, tracking_id_2) %>%
  pivot_longer(everything(), values_to = "tracking_id") %>% distinct(tracking_id)

# The grouping label is each feature's own putative name (or "m/z <mz>"), not a
# connected component of matched features, so two features linked only by mass
# stay separate rows when their annotations differ (e.g. D-lactate and threonic
# acid). The name is computed over the full significant set.
sig_tid <- sig %>% mutate(tracking_id = paste(compartment, timepoint, feature, sep = "__"))
group_lookup <- ct_name_track(sig_tid) %>% select(tracking_id, tracking_name)

d <- sig_tid %>%
  semi_join(tracked, by = "tracking_id") %>%
  left_join(group_lookup, by = "tracking_id") %>%
  mutate(compartment_timepoint = factor(paste(compartment, timepoint, sep = " · "),
                                        levels = CT_LEVELS),
         comp3     = unname(POOL3[compartment]),
         reliably  = is_true(supplement_reliably_detected))

# Inclusion is >= 2 distinct compartment x timepoint cells, as in the original
# analysis; a name can meet this within one compartment, so the across-compartment
# subset is counted separately below. Keep the lowest-P feature per name x cell.
up <- d %>% filter(direction == "Upregulated", !is.na(tracking_name), tracking_name != "") %>%
  group_by(tracking_name) %>%
  filter(n_distinct(compartment_timepoint) >= 2) %>%
  arrange(compartment_timepoint, raw_p, .by_group = TRUE) %>%
  distinct(tracking_name, compartment_timepoint, .keep_all = TRUE) %>%
  ungroup()

per <- up %>% group_by(tracking_name) %>%
  summarise(n_cells = n_distinct(compartment_timepoint),
            n_comp3 = n_distinct(comp3), reliably = any(reliably), .groups = "drop")
span2 <- per %>% filter(n_comp3 >= 2)

dir.create(dirname(ROWDATA), recursive = TRUE, showWarnings = FALSE)
write_csv(up, ROWDATA)
# unannotated "m/z" rows cannot be mapped to KEGG compounds
write_csv(tibble(compound = sort(per$tracking_name[!grepl("^m/z ", per$tracking_name) &
                                                   !per$tracking_name %in% JUNK_NAMES])),
          CMPDLIST)

cat("up-regulated names in >= 2 compartment x timepoint cells:", nrow(per), "\n")
cat("  spanning >= 2 of {maternal blood, milk, infant blood}:", nrow(span2), "\n")
cat("  ... of which reliably detected in the BEP supplement:", sum(span2$reliably),
    sprintf("(%.1f%%)\n", 100 * mean(span2$reliably)))
cat("  spanning all three:", sum(per$n_comp3 >= 3), "->",
    paste(sort(per$tracking_name[per$n_comp3 >= 3]), collapse = ", "), "\n")
cat("wrote", ROWDATA, "and", CMPDLIST, "\n")
