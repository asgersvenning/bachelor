library(MASS)
library(kernlab)
library(tidyverse)
library(magrittr)
library(hash)
library(asgerbachelor)
library(extrafont)

latin_df <- readxl::read_excel("../FIA_species_codes.xlsx") %>% 
  select(Genus, Species, `FIA Code`, `PLANTS Code`) %>% 
  mutate(
    name = paste(Genus,Species)
  )

fia_latin <- hash(latin_df$`FIA Code`, latin_df$name)
PLANTS_latin <- hash(latin_df$`PLANTS Code`, latin_df$name)
ECO_dict <- ecoregions_meta %>% 
  select(L2_KEY) %>% 
  distinct %>% 
  mutate(L2_KEY = str_split(L2_KEY,"  ") %>% 
           unlist %>% 
           matrix(nrow = 2) %>% 
           t %>% 
           set_colnames(c("key","value")) %>% 
           as_tibble) %>% 
  unnest(L2_KEY) %>% 
  {
    hash(.$key,.$value)
  }


eco <- ecoregionsToGrid(aggregationGrid, 2, 2500)

SPCoord <- FIA %>% 
  filter(SPCD %in% keys(fia_latin) & !is.na(ID) & ID %in% unique(eco2$ID)) %>% 
  mutate(name = values(fia_latin, SPCD)) %>% 
  group_by(x,y,ID,name) %>% 
  summarize(
    density = mean(individuals/samples),
    .groups = "drop"
  ) %>% 
  arrange(name) %>% 
  pivot_wider(c(x,y,ID),names_from = name, values_from = density, values_fill = 0)

eco_wide <- eco %>% 
  filter(ID %in% SPCoord$ID) %>%
  select(!geometry) %>% 
  arrange(2) %>% 
  pivot_wider(1,names_from=2,values_from=3,values_fill = 0) %>% 
  arrange(ID)

eco_classes <- eco_wide[,-1] %>% names %>% str_remove("eco") %>% factor

eco_max <- eco_wide[,-1] %>% 
  as.matrix %>% 
  apply(1, function(x) {
    x[is.na(x)] <- -Inf
    wm <- which(x == max(x)) 
    wm <- if (length(wm) == 1) wm else sample(wm, 1)
    eco_classes[wm]
  })

ecoCommunities <- SPCoord[,-c(1:2)] %>% 
  arrange(ID) %>% 
  mutate(across(!ID, ~.x/rowSums(across(!ID))),
         eco = factor(eco_max)) %>%
  # filter(between(rowSums(across(!c(ID,eco))>0),25,100)) %>% 
  relocate(eco,1)

ecoCommunitiesOrd <- ecoCommunities %>% 
  transmute(
    eco,
    cmdscale(dist(across(!c(eco,ID)),"manh")/2,18)
  ) %>% 
  rename("pcoa" = 2) %>% 
  mutate(pcoa = as_tibble(pcoa)) %>% 
  unnest(pcoa)


LDA_kfold <- lapply(1:5, function(rep) {
  folds <- ecoCommunitiesOrd %>% 
    group_by(eco) %>% 
    transmute(
      inTrain = sample(1:5,n(),T)
    ) %>% 
    pull(inTrain)
  
  lapply(1:5, function(x) {
    train <- x != folds
    
    LDA <- lda(eco ~ ., data = ecoCommunitiesOrd[train,])
    
    tibble(predicted = predict(LDA,ecoCommunitiesOrd[-train,])$class,
           truth = ecoCommunitiesOrd[-train,]$eco)
  }) %>% 
    {
      do.call(bind_rows,.)
    }
}) %>% 
  {
    do.call(bind_rows,.)
  }

LDA_kfold %>% 
  count(predicted,truth) %>% 
  complete(predicted,truth,fill=list(n=0)) %>% 
  group_by(truth) %>% 
  mutate(n = n/sum(n)) %>% 
  ungroup %>% 
  mutate(across(c(truth,predicted), ~hash::values(ECO_dict,
                                                  fct_reorder(.x, 
                                                              as.numeric(str_remove(.x,"\\.")))))) %>%
  ggplot(aes(predicted,truth,fill=n)) +
  geom_raster() +
  scale_fill_viridis_c(option = "A",
                       limits = 0:1,
                       labels = scales::percent) +
  labs(title = "Linear discriminant analysis",
       x = "Predicted class", y = "True class",
       fill = "Accuracy") +
  theme(text = element_text(family = "CMU Serif"),
        title = element_text(face = "bold",
                   %>%           size = 14),
        plot.title = element_text(hjust = .5,
                                  size = 20,
                                  face = "plain"),
        legend.key.width = unit(2,"lines"),
        legend.key.height = unit(1.5,"lines"),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        aspect.ratio = 1)


LDA_kfold %>% 
  count(predicted,truth) %>% 
  group_by(truth) %>% 
  mutate(accuracy = n/sum(n),
         w = n()) %>% 
  group_by(isTrue = truth == predicted,truth) %>% 
  summarize(accuracy = sum(accuracy),
            w = first(w)) %>% 
  # group_by(isTrue) %>% 
  # summarize(accuracy = weighted.mean(n,w)) %>% 
  filter(isTrue) %>% 
  pull(accuracy) %>% 
  round(3) %>% 
  as.character %>% 
  {
    write_lines(.,file("LDA_accuracy.txt"))
  }

