source("tidyData_and_calculate_indices.R")
library(stars)

baseTheme <-  ggpubr::theme_pubr() +
  theme(text = element_text(family = "CMU Serif"),
        title = element_text(face = "bold",
                             size = 14),
        strip.text = element_text(face = "bold",
                                  size = 14),
        legend.position = "right")

finalModelData <- modelDataRaw %>%  
  bind_cols(bind_rows(communityCoords,communityCoords)) %>% 
  mutate(Ecoregion = factor(Ecoregion),
         type = factor(type),
         Richness = log(Richness),
         across(c(Richness,`Shannon\nEvenness`),~scale(.x) %>% as.vector)) %>%
  rename("ShannonEvenness" = `Shannon\nEvenness`)


finalBiasModels <- sapply(c("FDis","FDiv","FEve"),function(x) gamm4::gamm4(formula(paste0(x, " ~ type + s(x,y,log(nPlots), by = type)")),
                                                                           random = ~ (1|Ecoregion/fid),
                                                                           data = finalModelData),
                          simplify = F,
                          USE.NAMES = T)


finalBiasResult <- tibble(results = finalBiasModels %>% 
                            lapply(function(x) summary(x$gam)$p.table %>% 
                                     as.data.frame %>% 
                                     rownames_to_column %>% 
                                     as_tibble)) %>% 
  mutate(Index = names(results)) %>% 
  unnest(results) %>% 
  filter(!str_detect(rowname,"Intercept")) %>% 
  select(!rowname)

write_csv(finalBiasResult,"finalBiasResult.csv")


