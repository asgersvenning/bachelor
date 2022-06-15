library(tidyverse)
library(magrittr)
library(hash)
library(asgerbachelor)
library(extrafont)

latin_df <- readxl::read_excel("FIA_species_codes.xlsx") %>% 
  select(Genus, Species, `FIA Code`, `PLANTS Code`) %>% 
  mutate(
    name = paste(Genus,Species)
  )

fia_latin <- hash(latin_df$`FIA Code`, latin_df$name)
PLANTS_latin <- hash(latin_df$`PLANTS Code`, latin_df$name)

SPCoord <- FIA %>% 
  filter(SPCD %in% keys(fia_latin) & !is.na(ID)) %>% 
  mutate(name = values(fia_latin, SPCD)) %>% 
  group_by(x,y,ID,name) %>% 
  summarize(
    density = mean(individuals/samples),
    .groups = "drop"
  ) %>% 
  arrange(name) %>% 
  pivot_wider(c(x,y,ID),names_from = name, values_from = density, values_fill = 0)

SPTraits <- PLANTS_traits %>% 
  mutate(name = values(PLANTS_latin, plants_code)) %>% 
  filter(name %in% names(SPCoord)[-c(1:3)]) %>% 
  select(!c(plants_code,
            contains("Propagated"),
            `Palatable Human`,
            `Protein Potential`,
            contains("Product"),
            `Planting Density per Acre, Maximum`,
            `Planting Density per Acre, Minimum`,
            `Hedge Tolerance`,
            `Commercial Availability`)) %>% 
  select(where(~mean(is.na(.x))<.33)) %>% 
  group_by(Genus = str_extract(name, "[:alpha:]+")) %>% 
  filter(mean(is.na(across(!name)))<.33) %>% 
  ungroup %>% 
  filter(rowMeans(is.na(across(!name)))<.33) %>% 
  relocate(name) %>% 
  arrange(name) %>%
  select(!Genus) %>% 
  gower_traits(T, NA.tolerance = 1)

SPTraits$traits <- SPTraits$traits %>% 
  group_by(Genus = str_extract(name, "[:alpha:]+")) %>% 
  mutate(across(!name, function(x) {
    replace(x,which(is.na(x)),mean(x,na.rm=T))
  })) %>% 
  ungroup %>% 
  select(!Genus)


SPCoord %>% 
  ggplot(aes(x,y)) +
  geom_point()

tBC <- SPCoord[,-1:-3] %>% 
  as.matrix %>% 
  divide_by(rowSums(.)) %>% 
  dist("manh") %>% 
  divide_by(2)

tCMD <- tBC %>% 
  cmdscale

cmdClust <- tBC %>% 
  hclust("ward.D2") %>% 
  cutree(7) %>% 
  factor

tCMD %>% 
  set_colnames(c("PC1","PC2")) %>% 
  as_tibble %>% 
  mutate(cluster = cmdClust) %>% 
  ggplot(aes(PC1,PC2, color = cluster)) +
  geom_point()

SPCoord %>% 
  # filter(rowSums(across(-c(1:3))!=0)>=10) %>% 
  mutate(cluster = cmdClust) %>% 
  ggplot(aes(x,y, fill = cluster)) +
  geom_raster() +
  scale_fill_brewer(palette = "Dark2") +
  theme_void() +
  coord_equal()


SPClust <- SPCoord %>% 
  mutate(cluster = cmdClust) %>% 
  relocate(cluster, .after = 3)


SPClust %>% 
  transmute(
    x,y,ID,cluster,
    rich = rowSums(across(!c(x,y,ID,cluster),~.x>0)),
    shEve = apply(across(!c(x,y,ID,cluster)), 1, 
                  function(x) {
                    -sum(x/sum(x) * log(x/sum(x)), na.rm=T)/log(sum(x>0))
                  })
  ) %>% 
  ggplot(aes(rich,shEve)) +
  geom_bin_2d(aes(color = after_stat(count))) +
  geom_smooth(method = "lm") +
  scale_fill_viridis_c(option = "A",
                       trans = "log10") +
  scale_color_viridis_c(option = "A",
                       trans = "log10",
                       guide = "none") +
  facet_wrap(~paste0("Cluster ", cluster), nrow = 2) +
  coord_cartesian(expand = F) +
  ggpubr::theme_pubr() +
  theme(
    aspect.ratio = 1,
    text = element_text(family = "CMU Serif"),
    title = element_text(face = "bold",
                         size = 14),
    panel.spacing.x = unit(2,"lines"),
    legend.position = "right"
  ) + 
  labs(
    y = "Shannon Evenness", x = "Species Richness"
  )
