tibble(models = finalBiasModels %>% lapply(function(x) x$gam)) %>% 
  mutate(residuals = map(models,function(x) tibble(fitted = fitted(x),
                                                   residual = residuals(x))),
         coords = bind_rows(communityCoords,communityCoords) %>% list,
         type = finalModelData %>%
           select(!c(x,y,L2_KEY)) %>% 
           list,
         Index = names(models)) %>%
  select(!models) %>% 
  unnest(everything()) %>% 
  nest(data = everything()) %>% 
  summarise(
    plots = c(fitted = map(data, ~ ggplot(.x, aes(fitted,residual, fill = type, color = type)) +
                             geom_bin_2d(aes(alpha = after_stat(count)),
                                         color = NA,
                                         key_glyph = draw_key_point) +
                             geom_hline(yintercept = 0, 
                                        color = "firebrick2", 
                                        linetype = "dashed", 
                                        size = 1) +
                             stat_summary_bin(fun.data = mean_sdl,
                                              bins = 20,
                                              show.legend = F) +
                             scale_color_brewer(palette = "Dark2",
                                                guide = "none") +
                             scale_fill_brewer(palette = "Dark2") +
                             scale_alpha_continuous(trans = "log10",
                                                    range = c(0.05,.75),
                                                    n.breaks = 6) +
                             facet_wrap(~Index, scales = "free") +
                             baseTheme +
                             guides(alpha = guide_legend(nrow = 2,
                                                         override.aes = list(shape = 21,
                                                                             size = 5)),
                                    fill = guide_legend(ncol = 1,
                                                        override.aes = list(shape = 21,
                                                                            size = 5))) +
                             theme(legend.box = "horizontal",
                                   legend.box.just = "right",
                                   legend.position = "top",
                                   aspect.ratio = 1) +
                             labs(x = "Fitted values", y = "Residual",
                                  alpha = "Density", fill = "Analysis\nMode")
    ),
    spatial = map(data, ~ .x %>%
                    select(x,y,residual,fitted,Index,type) %>% 
                    st_as_sf(coords = c("x","y")) %>% 
                    st_set_crs(st_crs(ecoregions_L2)) %>% 
                    ggplot(aes(fill=residual/abs(residual+fitted))) +
                    stat_sf_coordinates(aes(after_stat(x),after_stat(y)),
                                        geom = "tile") +
                    scale_fill_viridis_c(option = "A") +
                    facet_grid(rows = vars(Index), cols = vars(type)) +
                    coord_sf() +
                    baseTheme +
                    theme(panel.grid.major = element_line(color = "gray65",
                                                          size = .25,
                                                          linetype = "dashed")) +
                    labs(x = NULL, y = NULL, fill = latex2exp::TeX("$\\frac{\\hat{y}-y}{y}$"))
    ),
    ecoregion = map(data, ~ ggplot(.x,aes(residual,paste0(str_extract(Ecoregion,".{0,10}"),ifelse(str_detect(Ecoregion,".{11}"),"...","")))) +
                      geom_violin(scale = "width", size = 1) +
                      stat_summary(fun.data = mean_sdl,
                                   color = "red",
                                   geom = "errorbar",
                                   size = .75,
                                   width = .5) +
                      stat_summary(fun = mean,
                                   color = "red",
                                   fill = "gray65",
                                   shape = 21,
                                   geom = "point",
                                   stroke = 1,
                                   size = 4) +
                      geom_vline(xintercept = 0, color = "black", size = .75) +
                      baseTheme +
                      theme(axis.text.y = element_text(hjust = 0)) +
                      labs(y = NULL, x = "Residual") +
                      facet_wrap(~type,nrow=1,scales = "free_x")),
    samples = map(data, ~ ggplot(.x, aes(exp(nPlots),residual, fill = type, color = type)) +
                    geom_bin_2d(aes(alpha = after_stat(count)),
                                color = NA,
                                key_glyph = draw_key_point) +
                    geom_hline(yintercept = 0, 
                               color = "firebrick2", 
                               linetype = "dashed", 
                               size = 1) +
                    stat_summary_bin(fun.data = mean_sdl,
                                     bins = 20,
                                     show.legend = F) +
                    scale_color_brewer(palette = "Dark2",
                                       guide = "none") +
                    scale_fill_brewer(palette = "Dark2") +
                    scale_alpha_continuous(trans = "log10",
                                           range = c(0.05,.75),
                                           n.breaks = 6) +
                    scale_x_log10() +
                    facet_wrap(~Index, scales = "free") +
                    baseTheme +
                    guides(alpha = guide_legend(nrow = 2,
                                                override.aes = list(shape = 21,
                                                                    size = 5)),
                           fill = guide_legend(ncol = 1,
                                               override.aes = list(shape = 21,
                                                                   size = 5))) +
                    theme(legend.box = "horizontal",
                          legend.box.just = "right",
                          legend.position = "top",
                          aspect.ratio = 1) +
                    labs(x = "# Plots in Grid Cell", y = "Residual",
                         alpha = "Density", fill = "Analysis\nMode")
    )),
    names = names(plots),
    width = c(4,3.5,4,4),
    height = c(2,3,2.5,2)
  ) %>% 
  mutate(save = pmap(list(names,plots,width,height), function(type,plot,w,h) {
    ggsave(paste0(type,"_biasResidualPlot.png"),plot,
           type = "cairo-png", dpi = 400, scale = 3,
           width = w, height = h)
  }))


finalBiasResultPlot <- finalBiasResult %>% 
  rename(estimate = 1, std.error = 2, t.value = 3, p.value = 4) %>% 
  ggplot(aes(estimate,Index,xmin=estimate-std.error*1.96,xmax=estimate+std.error*1.96)) +
  geom_vline(xintercept = 0, color = "firebrick2", linetype = "dashed", size = 1) +
  geom_pointrange() +
  geom_mark_circle(aes(label = paste0("P-value = ", scales::scientific(p.value,2),
                                      "\nt-value  = ", str_extract(t.value,"[-]{0,1}.{5}"))),
                   # color = "#FFFFFF00",
                   expand = unit(3,"mm"),
                   radius = unit(3,"mm"),
                   label.buffer = unit(10,"mm"),
                   con.size = 1) +
  scale_x_continuous(expand = expansion(.35)) +
  scale_y_discrete(expand = expansion()) +
  coord_cartesian(ylim = c(.95,1.25)) +
  labs(y = NULL, x = latex2exp::TeX("$Index_{Abundance} - Index_{Presence/Absence}$")) +
  ggpubr::theme_pubr() +
  theme(text = element_text(family = "CMU Serif"),
        title = element_text(hjust = .5,
                             size = 20),
        strip.text = element_text(face = "bold",
                                  size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = .3) +
  facet_wrap(~Index, scales = "free", ncol = 1, strip.position = "bottom")

ggsave("finalBiasResult.png",finalBiasResultPlot,
       type = "cairo-png", dpi = 400, height = 2.4, width = 2, scale = 3)