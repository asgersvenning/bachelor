library(tidyverse)
library(ggforce)
library(ambient)
library(extrafont)


tibble(t = seq(0,2*pi,length.out=1000),
       x = cos(t),
       y = sin(t)) %>% 
  mutate(
    z = 2 + ambient::gen_perlin(x,y,frequency = 5,pertubation_amplitude=.5)) %>% 
  arrange(t) %>% 
  ggplot(aes(x*z,y*z)) +
  geom_point() 

rasterText <- function(string,cex=10,xlim=c(-1,1),ylim=c(-1,1),...) {
  temp <- tempfile(fileext = ".png")
  png(temp,type="cairo-png",...)
  plot.new()
  plot.window(xlim=xlim,ylim=ylim,mar=rep(0,4),oma=rep(0,4),xpd=NA)
  text(0,0,lab=string,cex=cex,family="CMU Serif")
  dev.off()

  png::readPNG(temp) %>% 
    apply(3,function(x) raster::raster(x)) %>% 
    raster::stack() %>% 
    sum %>% 
    raster::as.matrix() %>%
    apply(2,rev)
}


rasterTitle <- rasterText("Functional topology\n&\nNiche filling",20,height=500,width=1200,res=35) %>% 
  raster::raster() %>% 
  {
    . + terra::focal(.,matrix(1,ncol=5,nrow=5),sd,na.rm=T,expand=T,fillvalue=NA)  
  }

library(raster)
rasterTitle %>% 
  as.data.frame(xy=T) %>% 
  as_tibble %>% 
  filter(!is.na(layer) & layer < 2.5) %>%  
  ggplot(aes(x,-y,fill=gen_worley(x,y,frequency = .01,value="distance"))) +
  geom_tile(show.legend = F) +
  # stat_summary_hex(aes(color=after_scale(fill)),
  #                  bins = c(100,55),
  #                  show.legend = F,
  #                  fun = mean) +
  coord_cartesian(expand = T) +
  scale_fill_viridis_c(option = "A") +
  # scale_fill_gradient2(low = "white", mid = "gray95", high = "#AA1155") +
  theme_void() +
  theme(aspect.ratio = .45)


tibble(x = 1:250,
       y = x) %>% 
  complete(x,y) %>% 
  ggplot(aes(x,y,fill=gen_worley(x,y,frequency = 0.02,octave=5,value="distance"))) +
  geom_tile() +
  coord_cartesian(expand = T) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  theme(aspect.ratio = 1)
