library(rgdal)
library(grid)
library(extrafont)

.pt <- 1 / 0.352777778
len0_null <- function(x) {
  if (length(x) == 0)  NULL
  else                 x
}

theme_border <- function(
    type = c("left", "right", "bottom", "top", "none"),
    colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + opts( panel.border=theme_border(type=c("bottom","left")) ) + ...
  type <- match.arg(type, several.ok=TRUE)
  structure(
    list(type = type, colour = colour, size = size, linetype = linetype),
    class = c("theme_border", "element_blank", "element")
  )
}
element_grob.theme_border <- function(
    element, x = 0, y = 0, width = 1, height = 1,
    type = NULL,
    colour = NULL, size = NULL, linetype = NULL,
    ...) {
  if (is.null(type)) type = element$type
  xlist <- c()
  ylist <- c()
  idlist <- c()
  if ("bottom" %in% type) { # bottom
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y, y))
    idlist <- append(idlist, c(1,1))
  }
  if ("top" %in% type) { # top
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y+height, y+height))
    idlist <- append(idlist, c(2,2))
  }
  if ("left" %in% type) { # left
    xlist <- append(xlist, c(x, x))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(3,3))
  }
  if ("right" %in% type) { # right
    xlist <- append(xlist, c(x+width, x+width))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(4,4))
  }
  if (length(type)==0 || "none" %in% type) { # blank; cannot pass absence of coordinates, so pass a single point and use an invisible line
    xlist <- c(x,x)
    ylist <- c(y,y)
    idlist <- c(5,5)
    linetype <- "blank"
  }
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = xlist, y = ylist, id = idlist, ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}

cb <- as.vector(sapply(c("FDis","FDiv", "FEve", "FRic"),
                       function(x) paste0(c("Abundance","Presence/Absence"),x))) 

comGCoord <- communityCoords[,1:2] %>%
  as.matrix %>%
  rgdal::project(CRS(st_crs(ecoregions_L2)$wkt)@projargs, inv = T) %>%
  as_tibble

violinData <- modelDataRaw %>%
  bind_cols(bind_rows(comGCoord,comGCoord)) %>% 
  nest(dat = !Ecoregion) %>% 
  mutate(Ecoregion = str_extract(Ecoregion,"[[:alpha:]/-[ ]]+") %>% 
           str_replace("/"," / "),
         Ecoregion = map_chr(Ecoregion, str_wrap, width = 16),
         Ecoregion = map_chr(Ecoregion, function(x) {
           lines <- strsplit(x,"\n") %>% unlist
           if (length(lines)<=2) return(x)
           out <- lines[1:2]
           out[2] <- str_extract(out[2],".{1,12}")
           return(paste0(paste0(out,collapse="\n"),"..."))
         }))  %>% 
  unnest(dat) %>% 
  mutate(Ecoregion = fct_reorder(Ecoregion, x, .fun = mean)) %>% 
  arrange(Ecoregion) %>% 
  group_by(Ecoregion) %>% 
  mutate(x = mean(x)) %>%
  mutate(x = map_chr(abs(x), ~ifelse(.x < 0, paste(sprintf("%.0f", round(.x,0)), "°E"), ifelse(.x > 0, paste(sprintf("%.0f", round(.x,1)), "°W"),sprintf("%.0f", round(.x,0)))))) %>%
  ungroup %>% 
  pivot_longer(FDis:FRic,names_to="Index") %>% 
  group_by(Index) %>% 
  mutate(rline = ifelse(Index == "FRic", max(value), NA)) %>% 
  group_by(Index,Ecoregion) %>% 
  mutate(rline = ifelse(row_number() == 1, rline, NA))

indexEcoViolin <- violinData %>% 
  ggplot(aes(value,as.numeric(Ecoregion),fill=paste0(type,Index))) +
  geom_violin(scale = "width",
              key_glyph = draw_key_point,
              orientation = "y") +
  stat_summary(fun.data = mean_cl_normal,
               position = position_dodge(1),
               show.legend = F,
               fatten = 1.5,
               orientation = "y") +
  geom_vline(aes(xintercept = mean),
             data = modelDataRaw %>% 
               summarize(across(FDis:FRic,mean)) %>% 
               pivot_longer(FDis:FRic,names_to="Index",values_to="mean"),
             show.legend = F,
             color = "firebrick",
             linetype = "dashed",
             size = 1) +
  geom_vline(aes(xintercept = rline),
             size = 1) +
  scale_fill_brewer(palette = "Dark2",
                    labels = function(x) str_extract(x,"Abundance|Presence/Absence"),
                    breaks = cb) +
  scale_x_continuous(
    breaks = function(x) {
      ran <- range(x)
      rand <- abs(diff(ran))
      newran <- ran + rand*c(1,-1)*0.1
      signif(seq(newran[1],newran[2],length.out=3),2)
    },
    labels = function(x) str_remove_all(x,"0+(?=$)|(?<=^)0(?=\\.)")) +
  scale_y_reverse(breaks = violinData$Ecoregion %>% levels %>% length %>% seq,
  labels = violinData$Ecoregion %>% levels,
  sec.axis = dup_axis(labels = violinData %>% 
                        group_by(Ecoregion) %>% 
                        slice_head %>% 
                        pull(x))
  ) +
  facet_grid(cols = vars(Index), rows = vars(Ecoregion), scales = "free",switch="y") +
  coord_cartesian(expand = F) +
  guides(fill = guide_legend(override.aes = list(shape = 21,
                                                 size = 7,
                                                 color = "black",
                                                 stroke = 1),
                             nrow = 2,
                             keywidth = unit(2,"lines"))) +
  labs(fill = NULL) + 
  baseTheme +
  theme(axis.text.y = element_blank(),
        axis.text.y.right = element_text(family = "CMU Serif",
                                         hjust = 1),
        axis.line.y.right = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        strip.background.y = element_blank(),
        strip.background.x = element_rect(fill = "gray90",
                                          color = NA),
        strip.switch.pad.grid = unit(0,"lines"),
        strip.text.y.left = element_text(angle = 0,
                                         size = 9,
                                         hjust = 0),
        strip.text.x = element_text(size = 16),
        panel.border = theme_border("bottom"),
        panel.spacing.y = unit(0,"lines"),
        panel.spacing.x = unit(.5,"lines"),
        strip.placement = "outside",
        plot.margin = margin(.25,0,.25,.15,"lines"),
        legend.position = "bottom",
        legend.text = element_text(size = 11),
        legend.box.spacing = unit(0,"lines"))

indexEcoViolinGrob <- indexEcoViolin %>% 
  ggplotGrob

for(i in which(grepl("axis-r", indexEcoViolinGrob$layout$name))){
  indexEcoViolinGrob$grobs[[i]]$vp$x <- unit(0.85, "npc")     # originally 1npc
  indexEcoViolinGrob$grobs[[i]]$vp$valid.just <- c(1, 0.5) # originally c(1, 0.5
}

ggsave("index_by_ecoregion.png",indexEcoViolinGrob %>% gridExtra::grid.arrange(),
       type = "cairo-png", dpi = 400, width = 3.1, height = 2.8, scale = 3)


# tibble(x = 1:5,
#        lab = map_chr(x, ~paste0(letters[1:.x],collapse=""))) %>% 
#   slice_sample(n = 25, replace = T) %>% 
#   ggplot(aes(y=lab)) +
#   geom_bar() +
#   theme(axis.text.y = element_markdown(padding = unit(c(0,0,0,10), "mm"), debug = T),
#         aspect.ratio = .1,
#         strip.text = element_blank()) +
#   facet_wrap(~lab, scales = "free_y", ncol = 1)
