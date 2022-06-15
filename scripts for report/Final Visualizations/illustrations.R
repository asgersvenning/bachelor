# Counter-example

nPoints <- 4*10

noChangeFDis <- tibble(data = tibble(x = sin(seq(0,2*pi-2*pi/nPoints,2*pi/nPoints)),
                                     y = cos(seq(0,2*pi-2*pi/nPoints,2*pi/nPoints))) %>% 
                         list,
                       abundances = list(rep(1/nPoints,nPoints),(sin(seq(-4*pi,4*pi-pi*8/nPoints,pi*8/nPoints)) + 1.5) %>% 
                                           divide_by(sum(.))),
                       grp = c("Presence/Absence","Abundance")) %>% 
  unnest(everything()) %>%
  group_by(grp) %>% 
  mutate(
    centerX = round(weighted.mean(x,abundances),10),
    centerY = round(weighted.mean(y,abundances),10)
  ) %>% 
  ungroup %>% 
  ggplot(aes(x,y,color=abundances)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray65") +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "gray65") +
  geom_point(shape = 16,
             size = 6) +
  geom_point(inherit.aes = F,
             aes(centerX, centerY),
             color = "red",
             size = 10,
             shape = 16) + 
  geom_text_repel(inherit.aes = F,
                  aes(centerX,centerY,
                      label = "Centroid"),
                  force = 0,
                  max.overlaps = Inf,
                  color = "red",
                  size = 6,
                  fontface = "bold",
                  family = "CMU Serif",
                  nudge_y = .15,
                  nudge_x = .45) + 
  coord_equal() +
  scale_color_viridis_c(option = "A",
                        labels = c("Low","Mid","High"),
                        breaks = c(0.01, 0.025, 0.04)) + 
  facet_wrap(~grp, nrow = 1) +
  ggpubr::theme_pubr() +
  labs(x = "Trait 1",
       y = "Trait 2",
       color = "Relative\nWeight") +
  theme(legend.position = "right",
        text = element_text(family = "CMU Serif"),
        title = element_text(size = 16,
                             face = "bold"),
        strip.text = element_text(size = 16,
                                  face = "bold"),
        panel.spacing.x = unit(1,"lines"),
        plot.caption.position = "plot",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.width = unit(2.5,"lines"),
        legend.key.height = unit(2,"lines"))

ggsave("exampleNoChangeFDis.png",noChangeFDis,
       type = "cairo-png",width = 4, height = 2, scale = 2.35)

## East-West sampling imbalance
library(rgdal)



eastWestSampling <- modelData %>% 
  bind_cols(communityCoords[,1:2] %>%
              as.matrix %>%
              rgdal::project(CRS(st_crs(ecoregions_L2)$wkt)@projargs, inv = T) %>%
              as_tibble) %>%
  group_by(x_int = cut_interval(x,20,boundary = 0,dig.lab=10)) %>% 
  summarize(
    n = n(),
    p_mean = mean(nPlots),
    p = confint(glm(nPlots ~ 1,family = "quasipoisson"),trace = F) %>% 
      suppressMessages() %>% 
      exp %>% 
      set_names(c("min","max")) %>% 
      list
  ) %>% 
  mutate(x = map(x_int, function(x) {
    ran <- x %>% 
      str_extract_all("[[:digit:]\\.-]+") %>% 
      unlist %>% 
      as.numeric
    
    abs(c(min = ran[1], mid = mean(ran), max = ran[2]))
    
  })) %>% 
  unnest_wider(where(is.list),names_sep = "_") %>% 
  ggplot() +
  geom_rect(aes(xmin=x_min,xmax=x_max,ymin=0,ymax=n),
            fill = "#FFAAAA") +
  geom_pointrange(aes(x_mid,p_mean,ymax=p_max,ymin=p_min),
                  fatten = 3,
                  size = .75) +
  scale_y_log10(minor_breaks = seq(0,1000,50.01),
                breaks = c(10,30,50,100,300,500),
                sec.axis = dup_axis(name = "<p style='float: left;'># Grid Cells in Interval</p>")) +
  scale_x_reverse(labels = function(axis) sapply(axis, function(x) ifelse(x < 0, paste(x, "°E"), ifelse(x > 0, paste(x, "°W"),x)))) +
  coord_cartesian(expand = F,
                  ylim = c(10,700)) +
  baseTheme +
  theme(panel.grid.minor.y = element_line(color = "gray55",
                                          linetype = "dashed",
                                          size = .15),
        axis.title.y = element_markdown(),
        axis.title.y.right = element_markdown(),
        axis.line.y.right = element_line(colour = "#FF3333"),
        axis.ticks.y.right = element_line(colour = "#FF3333"),
        panel.background = element_rect(fill = NA),
        panel.ontop = T,
        aspect.ratio = 1,
        plot.caption = element_text(face = "plain")) +
  labs(x = "Latitude", y = "<p style='float: left;'># Plots in Grid Cell</p>")

ggsave("eastWestSampling.png",eastWestSampling,
       type = "cairo-png",dpi=400,width=4.5,height=4, scale = 1.5)

