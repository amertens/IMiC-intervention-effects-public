# Libraries

library(ggplot2)
library(readr)
library(tidyverse)

# Import data

adjusted_intervention_effects_results_clean <-
  readRDS(paste0(here::here(), "/results/adjusted_intervention_effects_results_clean.RDS"))

# Display data

head(adjusted_intervention_effects_results_clean)

# Generate heatmap

# Load required libraries
library(tidyverse)
library(ggplot2)



# Create the heatmap

# Generate the heatmap

adjusted_intervention_effects_results_clean %>%
  filter(str_starts(biomarker, "Tg."), measure == "ATE") %>%
  mutate(
    biomarker = factor(biomarker, levels = sort(unique(biomarker))),
    contrast = factor(contrast)
  ) %>%
  ggplot(aes(x = contrast, y = biomarker, fill = est)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#0571b0", mid = "white", high = "#ca0020", midpoint = 0,
    name = "ATE"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(title = "Average Treatment Effects (ATE) on Triglycerides")
