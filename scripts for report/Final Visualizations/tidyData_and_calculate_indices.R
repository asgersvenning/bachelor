library(asgerbachelor)
library(tidyverse)
library(magrittr)
library(ggrepel)
library(ggforce)
library(extrafont)
library(sf)

library(hash)

sapply(c("filter","select","summarise"),
       conflicted::conflict_prefer,
       winner = "dplyr") %>% 
  invisible
conflicted::conflict_prefer("values","hash")

# PLANTS_traits[,-1] %>%
#   summarize(across(everything(), function(x) {
#     vals <- unique(x[!is.na(x)])
#
#     out <- vals[1:min(c(5,length(vals)))]
#
#     if (length(out) < 5) c(out,rep(NA,5-length(out))) else out
#     })) %>%
#   t %>%
#   as.data.frame %>%
#   rownames_to_column %>%
#   as_tibble %>%
#   rename(variable = 1, value_1 = 2, value_2 = 3, value_3 = 4, value_4 = 5, value_5 = 6) %>%
#   mutate(type = NA) %>%
#   relocate(type, .after = 1) %>%
#   write_excel_csv("trait_types.csv")



trait_meta <- readxl::read_excel("trait_types.xlsx",
                   col_types = c("text","text","text","logical",
                                 "skip","skip","skip","skip","skip")) %>%
  mutate(type = map_chr(type, ~switch(.x,
                                      "F" = "Factor",
                                      "B" = "Binary",
                                      "C" = "Continous")),
         category = map_chr(category, ~switch(.x,
                                              "M" = "MORPHOLOGY/PHYSIOLOGY",
                                              "G" = "GROWTH REQUIREMENTS",
                                              "R" = "REPRODUCTION",
                                              "U" = "SUITABILITY/USE")),
         across(c(type,category),factor)) %>%
  mutate(n = map2_dbl(variable, type, function(x,y) {
    if (y != "Continous") length(unique(PLANTS_traits[[x]])) else NA
    }))

FIA_dict <- PLANTS_meta %>%
  arrange(plants_code) %>%
  filter(!str_detect(plants_code,"^[:digit:]")) %>%
  select(fia_code,plants_code) %>%
  {
    hash(.$plants_code,.$fia_code)
  }
PLANTS_dict <- invert(FIA_dict)

FIA_meta <- readxl::read_xlsx("C:/Users/asger/Desktop/6 Semester/Bachelor/Progress reports/Bachelor - Week 11/FIA_species_codes.xlsx") %>%
  select(`FIA Code`, Genus, Species,`PLANTS Code`) %>%
  mutate(latin = paste(Genus,Species)) %>%
  select(contains("Code"),latin) %>%
  filter(`PLANTS Code` %in% values(PLANTS_dict))

LATIN_dict <- hash(FIA_meta$`PLANTS Code`, FIA_meta$latin)
#
# # write_rds(LATIN_dict,"LATIN_dict.rds")
#
eco2_dict <- hash(ecoregions_meta$NA_L2CODE, ecoregions_meta$NA_L2NAME)

# eco2 <- ecoregionsToGrid(aggregationGrid, 2, 2500)
# 
# Clean and filter raw data
# 
# tSP <-
#   # read_csv2("FIA_grid10.csv") %>%
#   FIA %>%
#   filter(SPCD %in% values(FIA_dict)) %>%
#   select(!c(x,y)) %>%
#   nest(dat = !ID) %>%
#   mutate(dat = map(dat, function(x) {
#     x %>%
#       pivot_wider(c(INVYR,samples),names_from=SPCD,values_from=individuals,
#                   values_fill = 0) %>%
#       summarize(across(!c(INVYR,samples),~weighted.mean(.x,samples)),
#                 samples = sum(samples)) %>%
#       relocate(samples) %>%
#       pivot_longer(!samples,names_to="SPCD",values_to="individuals")
#   })) %>%
#   unnest(dat) %>%
#   mutate(plants_code = values(PLANTS_dict,SPCD)) %>%
#   select(!c(SPCD)) %>%
#   arrange(plants_code) %>%
#   pivot_wider(c(ID,samples),
#               names_from  = plants_code,
#               values_from = individuals,
#               values_fill = 0)
# 
# eco2_wide <- eco2 %>%
#   select(ID,L2_KEY,proportion) %>%
#   mutate(L2_KEY = as.character(L2_KEY) %>%
#            str_remove("[.]") %>%
#            factor) %>%
#   filter(ID %in% unique(tSP$ID)) %>%
#   mutate(eco = values(eco2_dict, L2_KEY)) %>%
#   select(!L2_KEY) %>%
#   arrange(eco) %>%
#   pivot_wider(ID,names_from=eco,values_from=proportion,
#               values_fill=0) %>%
#   arrange(ID)
# 
# tSP <- tSP %>%
#   filter(ID %in% eco2_wide$ID) %>%
#   arrange(ID) %>%
#   select(ID,samples,where(~any(.!=0)))
# 
# PLANTS_traits %>%
#   select(plants_code, all_of(trait_meta %>%
#                                filter(choose) %>%
#                                pull(variable))) %>%
#   filter(plants_code %in% names(tSP)[-c(1,2)]) %>%
#   select(!plants_code) %>%
#   transmute(nMissing = rowSums(across(everything(),is.na))) %>%
#   summarize(
#     mean = weighted.mean(nMissing,nMissing>0),
#     n = sum(nMissing > 0),
#     n_tot = n()
#   ) %>%
#   unlist %>%
#   paste0(collapse="\n") %>%
#   writeLines("traitCoverage.txt")
# 
# PLANTS_tG <- PLANTS_traits %>%
#   select(plants_code, all_of(trait_meta %>%
#                   filter(choose) %>%
#                   pull(variable) %>% 
#                     sort)) %>%
#   arrange(plants_code) %>%
#   gower_traits(T)
# 
# PLANTS_tG$traits <- PLANTS_tG$traits %>%
#   group_by(genus = values(LATIN_dict, plants.code) %>% str_extract("^[:alpha:]+(?= )")) %>%
#   mutate(across(!plants.code,function(x) {
#     x[is.na(x)] <- mean(x,na.rm=T)
#     x
#   })) %>%
#   ungroup %>%
#   select(!genus) %>%
#   filter(plants.code %in% names(tSP)[-1:-2])
# 
# PLANTS_tG$traits[,-1] %>% 
#   is.na %>% 
#   multiply_by_matrix(PLANTS_tG$weights) %>% 
#   as.vector %>% 
#   {
#     .[.!=0]
#   } %>% 
#   mean %>% 
#   as.character() %>% 
#   writeLines("afterGenusMean.txt")
# 
# PLANTS_tG$weights <- PLANTS_tG$weights[colMeans(is.na(PLANTS_tG$traits[,-1])) < 1/50]
# 
# PLANTS_tG$traits <- PLANTS_tG$traits %>%
#   select(where(~mean(is.na(.)) < 1/50)) %>%
#   filter(rowSums(across(!plants.code,~is.na(.x))) == 0)
# 
# writeLines(PLANTS_tG$weights %>% names %>% unique() %>% paste0(collapse=";"),
#            "usedTraits.txt")
# 
# tSP <- tSP %>%
#   select(ID,samples,all_of(PLANTS_tG$traits$plants.code)) %>%
#   filter(rowSums(across(!c(ID,samples),~.x>0))>=10)
# 
# communityCoords <- eco2 %>%
#   filter(ID %in% unique(tSP$ID)) %>%
#   group_by(ID) %>%
#   filter(proportion == max(proportion, na.rm=T)) %>%
#   ungroup %>%
#   arrange(ID) %>%
#   select(geometry,L2_KEY) %>%
#   unnest(geometry)
# 
# 
# TR <- PLANTS_tG$traits[,-1]
# TRW <- PLANTS_tG$weights
# fSP <- as.matrix(tSP[,-c(1:2)]/rowSums(tSP[,-c(1:2)]))
# 
# # Calculate functional indices
# 
# allIndices <- tibble(
#   communities = list(fSP,
#                      (fSP>0) %>%
#                        divide_by(rowSums(.))),
#   type = c("Abundance","Presence/Absence")
# ) %>%
#   mutate(
#     FDis = map(communities,~functional_dispersion(TR,TRW,a=.x)),
#     FDiv = map(communities,~functional_divergence(TR,TRW,a=.x,ch = 5)),
#     FEve = map(communities,~functional_evenness(TR,TRW,a=.x)),
#     FRic = map(communities,~functional_richness(TR,TRW,a=.x,ndim = 5))
#   ) %>%
#   select(!communities)
# 
# # Save objects to avoid having to recompute
# 
# saveRDS(allIndices,"allIndices.rds")
# saveRDS(TR,"TR.rds")
# saveRDS(TRW,"TRW.rds")
# saveRDS(tSP,"tSP.rds")
# saveRDS(fSP,"fSP.rds")
# saveRDS(communityCoords,"communityCoords.rds")
# 
# Load indices and tidied data from saved RDS files

allIndices <- readRDS("allIndices.rds")
TR <- readRDS("TR.rds")
TRW <- readRDS("TRW.rds")
tSP <- readRDS("tSP.rds")
fSP <- readRDS("fSP.rds")
communityCoords <- readRDS("communityCoords.rds")

# Combine tidied data and indicies, which will be used for modelling and visualizations

modelDataRaw <- allIndices %>% 
  mutate(Richness = rep(list(rowSums(fSP>0)),2),
         `Shannon\nEvenness` = rep(list(apply(fSP,1,function(x) -sum(x*log(x),na.rm=T)/log(sum(x>0)))),2),
         nPlots = rep(list(tSP$samples),2),
         Ecoregion = rep(list(values(eco2_dict,
                                     communityCoords$L2_KEY %>% 
                                       str_remove("[\\.]"))),2),
         fid = lapply(FDis,function(x) seq(length(x)))) %>% 
  unnest(everything()) %>% 
  mutate(fid = factor(fid)) 



modelData <- modelDataRaw %>% 
  summarise(
    across(c(FDis,FDiv,FEve,FRic),~.x %>% 
             matrix(ncol=2) %>% 
             apply(1, function(x) x[1] - x[2])
    ),
    across(!c(FDis,FDiv,FEve,FRic),~matrix(.x,ncol=2)[,1]
    )
  ) %>% 
  select(!type)
