library(tidyverse)
library(extrafont)

miamiVice <- c("#0BD3D3", "#F890E7")

indexConvergenceData <- modelDataRaw %>% 
  # mutate(FDis = fitted(finalREModels$FDis),
  #        FDiv = fitted(finalREModels$FDiv),
  #        FEve = fitted(finalREModels$FEve)) %>% 
  # select(!fid) %>% 
  # group_by(Ecoregion,type) %>% 
  # summarize(
  #   across(everything(),mean),
  #   n = n()
  # ) %>% 
  pivot_longer(FDis:FRic,
               names_to="Index",
               values_to="IndexValue") %>%
  pivot_longer(c(Richness,`Shannon\nEvenness`),
               names_to="CommunityStructure",
               values_to="Structure") %>% 
  mutate(
    Ecoregion = map_chr(Ecoregion, function(x) {
      if (nchar(x)>20) x <- paste0(str_extract(x,"^.{1,20}"),"...")
      x %>% str_wrap(15)
    })
  ) %>% 
  filter(Index != "FRic")


convergencePlots <- ggplot(indexConvergenceData, aes(Structure,IndexValue)) +
  geom_hex(aes(alpha=after_stat(count),
               fill = type,
               weight=nPlots),
           bins = 15) +
  geom_smooth(aes(group=type,weight=nPlots),
            color = "gray15",
            size = 1.75,
            show.legend = F,
            se = F) +
  geom_smooth(aes(color=type,weight=nPlots),
              size = .85,
              show.legend = F,
              se = F) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_alpha_continuous(range = c(0.05,.85),
                         n.breaks = 6) +
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(n.breaks = 6) +
  # scale_color_manual(values = miamiVice) +
  # scale_fill_manual(values = miamiVice) +
  facet_grid(rows = vars(Index), 
             cols = vars(CommunityStructure),
             scales="free",
             switch = "x") +
  coord_cartesian(expand = F) +
  guides(alpha = guide_legend(ncol = 2),
         fill = guide_legend(ncol = 1)) +
  labs(x = NULL,
       y = NULL,
       alpha = "Plots\nObserved",
       fill = "Analysis\nMode",
       color = NULL) +
  ggpubr::theme_pubr() +
  theme(text = element_text(family = "CMU Serif"),
        title = element_text(face = "bold",
                             size = 14),
        strip.text = element_text(face = "bold",
                                  size = 14),
        legend.position = "top",
        legend.box = "horizontal",
        legend.box.just = "right",
        aspect.ratio = 1,
        panel.spacing = unit(1.5,"lines")) 

ggsave("convergencePlot.png",convergencePlots,
       type = "cairo-png", dpi = 600, width = 4, height = 6, scale = 1.5)


indexDistributionPlot <- modelDataRaw %>% 
  select(FDis,FDiv,FEve,type) %>% 
  pivot_longer(!type,names_to="Index") %>% 
  ggplot(aes(Index,value,fill=type)) +
  geom_violin(size = .75,
              scale = "area",
              trim = F) +
  stat_summary(aes(group = paste0(Index,type)),
               fun.data = mean_sdl,
               fun.args = list(mult = 1),
               position = position_dodge(.9)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(expand = expansion(.5)) +
  scale_y_continuous(expand = expansion(.1),
                     breaks = seq(0,1,.1)) +
  labs(y = NULL, x = NULL, fill = NULL) +
  ggpubr::theme_pubr() +
  theme(text = element_text(family = "CMU Serif"),
        title = element_text(face = "bold",
                             size = 14),
        legend.text = element_text(face = "bold",
                                   size = 14),
        strip.text = element_text(face = "bold",
                                  size = 14),
        panel.spacing.x = unit(1.25,"lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 2) +
  facet_wrap(~Index, scales = "free", nrow = 1, strip.position = "bottom")

ggsave("indexDistribution.png",indexDistributionPlot,
       type = "cairo-png", dpi = 400, width = 4, height = 3, scale = 2)
