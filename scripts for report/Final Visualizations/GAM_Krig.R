library(mgcv)
library(sp)
library(stars)
data("meuse")

oS <- seq(100,750,50)
oN <- seq(10,75,5)
optP <- sapply(oN, function(n) sapply(oS, function(p) {
  bam(log10(zinc) ~ s(x,y,bs="gp",k=n,m=c(2,p)),
      data = meuse) %$%
    gcv.ubre
}))

tibble(
  GCV = optP %>%
    as.vector,
  n = rep(oN,each=length(oS)),
  p = rep(oS,length(oN))
) %>% 
  # mutate(GCV = GCV - min(GCV)) %>% 
  ggplot(aes(n,p,fill=GCV)) +
  geom_raster() +
  scale_fill_viridis_c(option = "A") +
  coord_cartesian(expand = F) +
  baseTheme +
  theme(aspect.ratio = 1)

tmod <- sapply(names(meuse)[3:6], function(x) bam(formula(paste0(x, " ~ s(x,y,bs='gp',k=45,m=c(2,500))")),
      data = meuse %>% 
        mutate(across(c(elev,dist),~as.vector(scale(.x)))),
      family = "quasipoisson"),USE.NAMES=T,simplify = F)

st_polygon(meuse %>% 
             select(x,y) %>% 
             as.matrix %>% 
             concaveman::concaveman(concavity = 2) %>% 
             list) %>% 
  st_sfc %>% 
  st_sf() %>% 
  st_buffer(50) %>% 
  st_rasterize(dx = 20, dy = 20) %>% 
  as_tibble %>% 
  drop_na %>% 
  select(!ID) %>% 
  mutate(elev = 1,
         dist = 1,
         soil = 1,
         lime = 1) %>%
  nest(dat = everything()) %>% 
  mutate(pred = lapply(tmod, mgcv::predict.bam, newdata = dat[[1]]) %>% 
           as_tibble %>% 
           list) %>% 
  unnest(everything()) %>% 
  mutate(across(cadmium:zinc,~as.vector(scale(.x)))) %>% 
  pivot_longer(cadmium:zinc,names_to="metal") %>% 
  ggplot(aes(x,y,fill=value)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_equal() +
  baseTheme +
  facet_wrap(~metal)

meuse %>% 
  ggplot(aes(dist.m,zinc)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()
