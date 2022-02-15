library(tidyverse)
library(sf)

## Load in data
# (The data is heavily preprocessed at this point!)
# The FIA data is aggregated two a grid, so no further processing is needed
FIA <- read_csv2(here::here("data-raw","FIA.csv"))
PLANTS <- read_csv2(here::here("data-raw","PLANTS.csv"))

# The PLANTS dataset is loaded in long format, and contains information
# which would be easier to handle as two separate tables:
# Connecting FIA codes with PLANTS codes
PLANTS_meta <- PLANTS %>% 
  select(fia_code, common_name, plants_code) %>% 
  distinct

# Species trait table
PLANTS_traits <- PLANTS %>% 
  select(plants_code, trait, value) %>% 
  group_by(plants_code, trait) %>% 
  summarize(
    value = first(value),
    .groups = "drop"
  ) %>% 
  pivot_wider(plants_code, 
              names_from = trait,
              values_from = value,
              values_fill = NA) %>% 
  select(!`NA`) 

# The EPA ecoregions are originally contained in a single collection of polygons
# which is disaggregated at the lowest resolution ("LEVEL 4"), but still annotated
# for all resolutions. This means that you will almost always want to dissolve some 
# of the polygons at different scales, and keeping the information about different
# resolutions is not important either (and it is also contained in the codes anyways,
# since the ecoregion coding scheme is hierarchichal). I have thus predissolved and/or
# removed unnessary attributes for 5 different levels.

# Level 0 -  NO regions
ecoregions_L0 <- read_sf(here::here("data-raw","EPA","ecoregions_boundary.shp"))

# Level 1
ecoregions_L1 <- read_sf(here::here("data-raw","EPA","dissolved","ecoregions_L1.shp"))

# Level 2
ecoregions_L2 <- read_sf(here::here("data-raw","EPA","dissolved","ecoregions_L2.shp"))

# Level 3
ecoregions_L3 <- read_sf(here::here("data-raw","EPA","dissolved","ecoregions_L3.shp"))

# Level 4
ecoregions_L4 <- read_sf(here::here("data-raw","EPA","dissolved","ecoregions_L4.shp"))

## Saving the data
# FIA
usethis::use_data(FIA, overwrite = T)

# PLANTS meta
usethis::use_data(PLANTS_meta, overwrite = T)

# PLANTS traits
usethis::use_data(PLANTS_traits, overwrite = T)

# Ecoregions level 0
usethis::use_data(ecoregions_L0, overwrite = T)

# Ecoregions level 1
usethis::use_data(ecoregions_L1, overwrite = T)

# Ecoregions level 2
usethis::use_data(ecoregions_L2, overwrite = T)

# Ecoregions level 3
usethis::use_data(ecoregions_L3, overwrite = T)

# Ecoregions level 4
usethis::use_data(ecoregions_L4, overwrite = T)
