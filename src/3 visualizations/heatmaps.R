
#use code from hackathon


rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
library(superheat)

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

#c(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit)

plotdf <- d  %>% filter(study=="Misame") %>% select("arm", all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit)
plotdf <- predict(preProcess(plotdf, method = c("medianImpute", "nzv")),plotdf) %>% as.data.frame()

#Heatmap
p1 <- superheat(plotdf[,-1], 
                scale=TRUE,
                col.dendrogram = TRUE,
                row.dendrogram = FALSE,
                membership.rows = plotdf$arm, #update to classify treatment arms
                title="",
                pretty.order.rows = TRUE,
                pretty.order.cols = FALSE)
plot(p1$plot)


plotdf2 <- scale(plotdf[,-1])
p1 <- superheat(plotdf2, 
                scale=F)
plot(p1$plot)


plotdf <- d  %>% filter(study=="Misame") %>% select("arm", all_milk_components$metabolomics[grepl("tg.",all_milk_components$metabolomics)])
plotdf <- predict(preProcess(plotdf, method = c("medianImpute", "nzv")),plotdf) %>% as.data.frame()
plotdf2 <- scale(plotdf[,-1])

#Heatmap
p3 <- superheat(plotdf2, 
                scale=TRUE,
                col.dendrogram = FALSE,
                row.dendrogram = FALSE,
                membership.rows = plotdf$arm, #update to classify treatment arms
                title="",
                pretty.order.rows = TRUE,
                pretty.order.cols = TRUE)
plot(p3$plot)



saveRDS(list(heatmap_misame_primary=p1$plot), file=paste0(here::here(),"/figures/figure-data/heatmaps.RDS"))
