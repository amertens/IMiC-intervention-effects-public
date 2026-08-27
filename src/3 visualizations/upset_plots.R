# Specific library
library(UpSetR)

resfull <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))

res <- resfull %>% filter(sigFDR==1)
unique(res$studytime)

res_list <- list(
  `Elicit (1 mo.)` =res$biomarker[res$studytime=="Elicit (1 mo.)"],
  `Elicit (5 mo.)` =res$biomarker[res$studytime=="Elicit (5 mo.)"],
  `Misame (14-21 days)` =res$biomarker[res$studytime=="Misame (14-21 days)"],
  `Misame (1-2 mo.)` =res$biomarker[res$studytime=="Misame (1-2 mo.)"],
  `Misame (3-4 mo.)` =res$biomarker[res$studytime=="Misame (3-4 mo.)"],
  `Vital (1.5 mo.)` =res$biomarker[res$studytime=="Vital (1.5 mo.)"],
  `Vital (2 mo.)` =res$biomarker[res$studytime=="Vital (2 mo.)"])
input = fromList(res_list)
set_order <- rev(names(res_list))

# Find elements that appear in all sets
common_elements <- Reduce(intersect, res_list)

# Plot
full_plot <- upset(input, 
      nintersects = 40, 
      nsets = 7,  # Matches the number of sets in res_list
      sets = set_order,  # This ensures the sets are ordered as they appear in res_list
      keep.order = TRUE, # This keeps the order specified in 'sets'
      order.by = "freq", # Orders intersections by frequency
      decreasing = TRUE, # Ensures most frequent intersections are on the left
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1) 


res_list <- list(
  `Elicit (1 mo.)` =res$biomarker[res$studytime=="Elicit (1 mo.)" & res$outcome_group != "tertiary"],
  `Elicit (5 mo.)` =res$biomarker[res$studytime=="Elicit (5 mo.)" & res$outcome_group != "tertiary"],
  `Misame (14-21 days)` =res$biomarker[res$studytime=="Misame (14-21 days)" & res$outcome_group != "tertiary"],
  `Misame (1-2 mo.)` =res$biomarker[res$studytime=="Misame (1-2 mo.)" & res$outcome_group != "tertiary"],
  `Misame (3-4 mo.)` =res$biomarker[res$studytime=="Misame (3-4 mo.)" & res$outcome_group != "tertiary"],
  `Vital (1.5 mo.)` =res$biomarker[res$studytime=="Vital (1.5 mo.)" & res$outcome_group != "tertiary"],
  `Vital (2 mo.)` =res$biomarker[res$studytime=="Vital (2 mo.)" & res$outcome_group != "tertiary"])
input = fromList(res_list[3:7])
set_order <- rev(names(res_list[3:7]))

# Find elements that appear in all sets
common_elements <- Reduce(intersect, res_list[3:5])
res_labs= res %>% distinct(category, biomarker)
res_labs$category[res_labs$biomarker %in% common_elements]

common_elements <- Reduce(intersect, res_list[3:7])
res_labs= res %>% distinct(category, biomarker)
res_labs$category[res_labs$biomarker %in% common_elements]

# Plot
misame_vital_plot <- upset(input, 
      nintersects = 25, 
      nsets = 5,  # Matches the number of sets in res_list
      sets = set_order,  # This ensures the sets are ordered as they appear in res_list
      keep.order = TRUE, # This keeps the order specified in 'sets'
      order.by = "freq", # Orders intersections by frequency
      decreasing = TRUE, # Ensures most frequent intersections are on the left
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1) 

#ggsave doesnt work
# ggsave(paste0(here::here(),"/figures/upset_plot_full.png"), plot=full_plot, width=10, height=6)
# ggsave(paste0(here::here(),"/figures/upset_plot_misame_vital.png"), plot=misame_vital_plot, width=10, height=6)


#how to add scatter plots below:

# movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )
# upset(movies,attribute.plots=list(gridrows=60,plots=list(list(plot=scatter_plot, x="ReleaseDate", y="AvgRating"),
#                                                          list(plot=scatter_plot, x="ReleaseDate", y="Watches"),list(plot=scatter_plot, x="Watches", y="AvgRating"),
#                                                          list(plot=histogram, x="ReleaseDate")), ncols = 2))


