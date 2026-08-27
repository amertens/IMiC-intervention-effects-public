
# Facet Plot

## Synthesize data


# Load necessary libraries
library(tidyverse)

# Set seed for reproducibility
set.seed(123)

# Function to generate synthetic data for a given time point, study, assay, and p-value
generateData <- function(timePointLabel, studyLabel, assayLabel, interventionMean, controlMean, pValue) {
  tibble(
    id = 1:30,
    concentration = rnorm(30, mean = interventionMean, sd = 100),  # Intervention group values
    group = "intervention",
    timePoint = timePointLabel,
    study = studyLabel,    # Assign study label
    assay = assayLabel,    # Assign assay label
    pValue = pValue        # Assign p-value
  ) %>%
    bind_rows(
      tibble(
        id = 31:60,
        concentration = rnorm(30, mean = controlMean, sd = 100),  # Control group values
        group = "control",
        timePoint = timePointLabel,
        study = studyLabel,    # Assign study label
        assay = assayLabel,    # Assign assay label
        pValue = pValue        # Assign p-value
      )
    )
}

# Define p-values for each time point (ensuring all are < 0.05)
pValues <- c("14-21" = 0.032, "1-2" = 0.0018, "3-4" = 0.00045, "1" = 0.029, "5" = 0.041, "1.5" = 0.037, "2" = 0.033)

# Generate datasets for each study, time point, and assay
syntheticData <- bind_rows(
  # MISAME-III Study (Original) for b1
  generateData("14-21", "misame3", "b1", 300, 200, pValues["14-21"]),
  generateData("1-2", "misame3", "b1", 300, 200, pValues["1-2"]),
  generateData("3-4", "misame3", "b1", 300, 200, pValues["3-4"]),
  
  # elicit Study for b1
  generateData("1", "elicit", "b1", 300, 500, pValues["1"]),
  generateData("5", "elicit", "b1", 300, 500, pValues["5"]),
  
  # mumpta Study for b1
  generateData("1.5", "mumpta", "b1", 300, 250, pValues["1.5"]),
  generateData("2", "mumpta", "b1", 300, 250, pValues["2"]),
  
  # MISAME-III Study for b2
  generateData("14-21", "misame3", "b2", 600, 500, pValues["14-21"]),
  generateData("1-2", "misame3", "b2", 600, 500, pValues["1-2"]),
  generateData("3-4", "misame3", "b2", 600, 500, pValues["3-4"]),
  
  # elicit Study for b2
  generateData("1", "elicit", "b2", 700, 600, pValues["1"]),
  generateData("5", "elicit", "b2", 300, 200, pValues["5"]),
  
  # mumpta Study for b2
  generateData("1.5", "mumpta", "b2", 600, 400, pValues["1.5"]),
  generateData("2", "mumpta", "b2", 600, 400, pValues["2"]),
  
  # MISAME-III Study for b3
  generateData("14-21", "misame3", "b3", 500, 150, pValues["14-21"]),
  generateData("1-2", "misame3", "b3", 500, 150, pValues["1-2"]),
  generateData("3-4", "misame3", "b3", 500, 150, pValues["3-4"]),
  
  # elicit Study for b3
  generateData("1", "elicit", "b3", 500, 150, pValues["1"]),
  generateData("5", "elicit", "b3", 500, 150, pValues["5"]),
  
  # mumpta Study for b3
  generateData("1.5", "mumpta", "b3", 500, 150, pValues["1.5"]),
  generateData("2", "mumpta", "b3", 500, 150, pValues["2"]),
  
  # MISAME-III Study for b6
  generateData("14-21", "misame3", "b6", 900, 700, pValues["14-21"]),
  generateData("1-2", "misame3", "b6", 900, 700, pValues["1-2"]),
  generateData("3-4", "misame3", "b6", 900, 700, pValues["3-4"]),
  
  # elicit Study for b6
  generateData("1", "elicit", "b6", 700, 700, pValues["1"]),
  generateData("5", "elicit", "b6", 800, 700, pValues["5"]),
  
  # mumpta Study for b6
  generateData("1.5", "mumpta", "b6", 800, 700, pValues["1.5"]),
  generateData("2", "mumpta", "b6", 1000, 700, pValues["2"])
)

# Print first few rows
print(syntheticData)

# Save dataset as CSV
write_csv(syntheticData, "syntheticData_extended.csv")

# Define a function to convert p-values to significance asterisks
getAsterisks <- function(p) {
  if(p < 0.001) {
    return("***")
  } else if(p < 0.01) {
    return("**")
  } else if(p < 0.05) {
    return("*")
  } else {
    return("")
  }
}

# Apply the function to each row
syntheticData <- syntheticData %>%
  mutate(asterisks = map_chr(pValue, getAsterisks))

### STEP 2: Your faceted box plot code (with a few adjustments)

# Define cutoff values for a subset of assays (here: B1, B2, B3, B6)
cutoffValues <- tibble(
  assay = c("B1", "B2", "B3", "B6"),  
  cutoff = c(400, 500, 400, 800)  
)

# Compute the maximum concentration per assay and set the y-position for asterisks
asteriskData <- syntheticData %>%
  group_by(assay) %>%
  summarise(
    maxConcentration = max(concentration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(asteriskY = maxConcentration * 1.07)

# Merge the asterisk data back to include studyTimeLabel (for positioning in each facet)
asteriskData <- syntheticData %>%
  distinct(studyTimeLabel, assay, asterisks) %>%
  left_join(asteriskData, by = "assay")

# Generate random percentage labels for each group within each boxplot.
# (Here, we compute the median concentration per studyTimeLabel, assay, and group,
# then sample a random percentage between 10 and 99.)
set.seed(123)
percentLabels <- syntheticData %>%
  group_by(studyTimeLabel, assay, group) %>%
  summarise(
    medianY = median(concentration, na.rm = TRUE),
    randomPercent = sample(10:99, 1),
    .groups = "drop"
  )

# Create the faceted box plot
facetedBoxPlot <- 
  ggplot(syntheticData, aes(x = studyTimeLabel, 
                            y = concentration, 
                            fill = group)) +
  # Add horizontal cutoff lines (only for assays that match cutoffValues)
  geom_hline(  
    data = cutoffValues,  
    aes(yintercept = cutoff),  
    color = "green", linetype = "solid", linewidth = 1  
  ) +
  # Add jittered points
  geom_jitter(
    aes(color = group),  
    position = position_jitterdodge(),  
    size = 1.5, alpha = 0.6  
  ) +  
  # Add boxplots (with outliers removed)
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  
  # Facet the plot by assay; each assay gets its own y-scale
  facet_wrap(~ assay, scales = "free_y") +  
  scale_fill_manual(values = c("intervention" = "red", "control" = "blue")) +  
  scale_color_manual(values = c("intervention" = "red", "control" = "blue")) +  
  labs(
    title = "Faceted Box Plot of Concentrations by Study and Time Point",
    x = "Study · Time Point",
    y = "Concentration",
    fill = "Group"
  )

# Add dynamically positioned asterisks (using the precomputed positions)
facetedBoxPlot <- facetedBoxPlot + 
  geom_text(
    data = asteriskData,  
    aes(x = studyTimeLabel, 
        y = asteriskY,  
        label = asterisks),  
    inherit.aes = FALSE,  
    size = 5
  )

# Add percentage labels (using geom_label) positioned at the median for each group
facetedBoxPlot <- facetedBoxPlot + 
  geom_label(
    data = percentLabels,  
    aes(x = studyTimeLabel, 
        y = medianY,  
        label = randomPercent),  
    inherit.aes = FALSE,  
    size = 2.5, fontface = "bold", 
    color = "white", fill = "black", label.size = 0.3, 
    label.r = unit(0.15, "lines")
  )

# Print the plot
print(facetedBoxPlot)

# # Save the plot as an A4 square image (210 x 210 mm)
# ggsave(
#   filename = "../results/facetedBoxPlot_A4_Square.png",  
#   plot = facetedBoxPlot,
#   device = "png",
#   width = 400,   # width in mm
#   height = 2000,  # height in mm
#   limitsize = FALSE,
#   units = "mm",
#   dpi = 600      # high resolution
# )





























































