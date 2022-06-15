source("tidyData_and_calculate_indices.R")

library(tidyverse)
library(magrittr)
library(ggrepel)
library(ggforce)
library(extrafont)
library(lme4)
library(patchwork)
# library(spatstat)
library(sf)
library(stars)
library(terra)
library(raster)
library(concaveman)

modelDataAdjusted <- modelDataRaw


visCRS <- st_crs(ecoregions_L2)

outlineSMPL <- ecoregions_L0 %>%
  st_transform(visCRS) %>% 
  st_simplify(dTolerance = 1000) %>% 
  st_buffer(2500) %>% 
  st_simplify(dTolerance = 1000)

wdim <- c(10,10)
mult <- 3
weights <- sapply(seq(-mult,mult,length.out=wdim[1]*2+1), 
                  function(x) sapply(seq(-mult,mult,length.out=wdim[2]*2+1), 
                                     function(y) 1/(1+sqrt(x^2+y^2)))) %>% 
  matrix(wdim[2]*2+1, wdim[1]*2+1) %>% 
  replace(which(.<quantile(.,.35)), NA)

wmid <- which(weights == 1)

## The next section has to be run in order from start to finish

modelDataInterp <- modelDataRaw %>%
  group_by(type) %>% 
  mutate(coord = communityCoords) %>% 
  unnest(coord)  %>% 
  ungroup %>% 
  pivot_longer(FDis:FRic,names_to="Index",values_to="Index_value") %>%
  group_by(type,Index) %>% 
  arrange(x,y) %>% 
  ungroup %>% 
  nest(dat = !c(type,Index)) %>% 
  mutate(
    dat = map(dat, function(xdat) {
      xdat %>% 
        dplyr::select(x,y,Index_value) %>% 
        mutate(across(c(x,y), ~.x - 25000)) %>%
        st_as_sf(coords = c("x","y")) %>% 
        st_rasterize(dx = 50000, dy = 50000, pretty=T,inside=F,
                     xlim = st_bbox(outlineSMPL)[c("xmin", "xmax")],
                     ylim = st_bbox(outlineSMPL)[c("ymin", "ymax")]) %>% 
        as("Raster")  %>% 
        extend(wdim+10) %>%
        terra::focal(
          w=matrix(1,wdim[2]*2+1,wdim[1]*2+1),
          fun=function(x) {
            
            mid <- x[wmid]
            
            if (!is.na(mid)) return(mid)
            if (all(is.na(x))) return(NA)
            
            w <- weights[!is.na(weights) & !is.na(x)]
            w <- w/sum(w)
            x <- x[!is.na(x) & !is.na(weights)]
            
            sum(x * w, na.rm=T)
          }) %>% 
        as("SpatRaster") %>% 
        as.data.frame(xy = T, na.rm=T) %>% 
        as_tibble %>% 
        rename(Index_value = 3)
    })
  ) %>% 
  mutate(
    Index = Index %>% 
      replace(which(Index == "FDis"),
              "Functional Dispersion") %>% 
      replace(which(Index == "FDiv"),
              "Functional Divergence") %>% 
      replace(which(Index == "FEve"),
              "Functional Evenness") %>% 
      replace(which(Index == "FRic"),
              "Functional Richness") %>% 
      factor
  ) 

newAgr <- modelDataInterp$dat[[1]] %>%
  select(x,y) %>% 
  mutate(fid = row_number()) %>% 
  st_as_stars(coords = c("x","y")) %>%
  st_as_sf %>%  
  st_intersection(outlineSMPL %>% 
                    st_set_crs(NA)) %>% 
  st_set_crs(visCRS) %>% 
  arrange(fid)

modelDataInterp <- modelDataInterp %>% 
  mutate(dat = map(dat, function(x) {
    newAgr %>% 
      select() %>% 
      mutate(value = x %>% 
               filter(row_number() %in% newAgr$fid) %>% 
               pull(Index_value))
  }))

modelDataNest <- modelDataRaw %>%
  group_by(type) %>% 
  mutate(coord = communityCoords) %>% 
  unnest(coord)  %>% 
  ungroup %>% 
  pivot_longer(FDis:FRic,names_to="Index",values_to="Index_value") %>%
  group_by(type,Index) %>% 
  arrange(x,y) %>% 
  ungroup %>% 
  nest(dat = !c(type,Index)) %>% 
  mutate(
    Index = Index %>% 
      replace(which(Index == "FDis"),
              "Functional Dispersion") %>% 
      replace(which(Index == "FDiv"),
              "Functional Divergence") %>% 
      replace(which(Index == "FEve"),
              "Functional Evenness") %>% 
      replace(which(Index == "FRic"),
              "Functional Richness") %>% 
      factor
  ) 

oldAgr <- modelDataNest$dat[[1]] %>%
  select(x,y) %>%  
  mutate(fid = row_number()) %>%
  mutate(across(c(x,y), ~.x - 25000)) %>%
  st_as_sf(coords = c("x","y")) %>% 
  st_rasterize(dx = 50000, dy = 50000, pretty=T,inside=F,
               xlim = st_bbox(outlineSMPL)[c("xmin", "xmax")],
               ylim = st_bbox(outlineSMPL)[c("ymin", "ymax")]) %>%
  st_as_sf %>%  
  st_intersection(outlineSMPL %>% 
                    select() %>% 
                    st_set_crs(NA)) %>%
  arrange(fid) %>% 
  st_set_crs(visCRS)

modelDataNest <- modelDataNest %>% 
  mutate(dat = map(dat, function(x) {
    oldAgr %>% 
      select() %>% 
      mutate(value = x %>% 
               arrange(x,y) %>% 
               filter(row_number() %in% oldAgr$fid) %>% 
               pull(Index_value))
  }))

## Plotting

plot_map <- function(x,y,z,lab,lab_x,lab_y) {
  IndexRange <- modelDataInterp %>%
    filter(Index == y) %>% 
    mutate(value = map(dat, ~.x$value)) %>% 
    pull(value) %>% 
    range
  
  ggplot(z, aes(fill=value)) +
    geom_sf(color = NA) +
    geom_sf(data = outlineSMPL,
            color = "black",
            fill = "#FFFFFF00",
            size = .25,
            inherit.aes = F) +
    annotate("text",x=lab_x,y=lab_y,label=lab,
             family = "CMU Serif",
             fontface = "bold",
             size = 6) +
    coord_sf() +
    scale_fill_viridis_c(option = c("A","D","E","C")[y],
                         limits = IndexRange) +
    labs(title = if (str_detect(y,"Dispersion")) x else NULL,
         fill = y %>% str_replace(" ", "\n"),
         x = NULL,
         y = NULL) +
    ggpubr::theme_pubr() + 
    theme(text = element_text(family = "CMU Serif"),
          plot.title = element_text(size = 24,
                                    hjust = .5,
                                    face = "bold"),
          legend.title = element_text(face = "bold",
                                      size = 14),
          legend.position = "right",
          plot.margin = margin(),
          axis.line.y = if (x == "Abundance") element_line() else element_blank(),
          axis.ticks.y = if (x == "Abundance") element_line() else element_blank(),
          axis.text.y = if (x == "Abundance") element_text() else element_blank(),
          axis.line.x = if (y == "Functional Richness") element_line() else element_blank(),
          axis.ticks.x = if (y == "Functional Richness") element_line() else element_blank(),
          axis.text.x = if (y == "Functional Richness") element_text() else element_blank()
    )
}

interpolatedMaps <- modelDataInterp %>% 
  mutate(annLab = paste0("(",letters[Index],".",ifelse(type == "Abundance",1,2),")"),
         ann_x = st_bbox(outlineSMPL)["xmin"] + 0.05*diff(st_bbox(outlineSMPL)[c("xmin","xmax")]),
         ann_y = st_bbox(outlineSMPL)["ymin"] + 0.05*diff(st_bbox(outlineSMPL)[c("ymin","ymax")])) %>% 
  mutate(maps = pmap(list(type, Index, dat, annLab, ann_x, ann_y), plot_map)) %>% 
  group_by(Index) %>%
  summarize(
    maps = wrap_plots(maps,guides="collect") %>% list
  ) %>% 
  pull(maps) %>% 
  wrap_plots(ncol = 1) 
# plot_annotation(title = "Data type",
#                 theme = theme(plot.title = element_text(family = "CMU Serif",
#                                                         face = "bold",
#                                                         size = 24,
#                                                         hjust = .5)))

ggsave("interpolatedMap.png",interpolatedMaps,
       type = "cairo-png", dpi = 600, width = 4, height = 4.5, scale = 3)

realMap <- modelDataNest %>% 
  mutate(annLab = paste0("(",letters[Index],".",ifelse(type == "Abundance",1,2),")"),
         ann_x = st_bbox(outlineSMPL)["xmin"] + 0.05*diff(st_bbox(outlineSMPL)[c("xmin","xmax")]),
         ann_y = st_bbox(outlineSMPL)["ymin"] + 0.05*diff(st_bbox(outlineSMPL)[c("ymin","ymax")])) %>% 
  mutate(maps = pmap(list(type, Index, dat, annLab, ann_x, ann_y), plot_map)) %>% 
  group_by(Index) %>%
  summarize(
    maps = wrap_plots(maps,guides="collect") %>% list
  ) %>% 
  pull(maps) %>% 
  wrap_plots(ncol = 1) 
# plot_annotation(title = "Data type",
#                 theme = theme(plot.title = element_text(family = "CMU Serif",
#                                                         face = "bold",
#                                                         size = 24,
#                                                         hjust = .5)))

ggsave("realMap.png", realMap,
       type = "cairo-png", dpi = 600, width = 4, height = 4.5, scale = 3)
