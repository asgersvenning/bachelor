library(spdep)

lw <- bind_rows(communityCoords,communityCoords)[,1:2] %>%
  st_as_sf(coords = 1:2) %>% 
  st_set_crs(st_crs(ecoregions_L2)) %>% 
  dnearneigh(0,5*10^5) %>% 
  nb2listw(zero.policy = T)



moranTests <- lapply(finalREModels, function(x) moran.test(residuals(x$gam),lw))

library(ncf)

splineDat <- communityCoords %>% 
  mutate(v = residuals(finalREModels$FDis$gam) %>% 
           # {
           #   .[(length(.)/2 + 1):length(.)]
           # }
           {
             .[1:(length(.)/2)]
           }
         ) 

spline.correlog(splineDat$x,splineDat$y,splineDat$v,
                resamp = 20,
                npoints = 100,
                xmax = 250000)
