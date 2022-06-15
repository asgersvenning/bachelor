tibble(models = finalREModels %>% lapply(function(x) x$mer)) %>% 
  mutate(residuals = map(models,function(x) tibble(fitted = fitted(x),
                                                   residual = residuals(x))),
         coords = bind_rows(communityCoords,communityCoords) %>% list,
         type = finalModelData %>%
           select(!c(x,y,L2_KEY)) %>% 
           list,
         Index = names(models)) %>%
  select(!models) %>% 
  unnest(everything()) %>% 
  select(residual,fid,type,Index) %>% 
  pivot_wider(c(fid,Index),names_from=type,values_from=residual) %>% 
  ggplot(aes(Abundance,`Presence/Absence`,fill=Index)) +
  geom_bin_2d(aes(alpha = after_stat(count))) +
  geom_abline(slope = 1, size = 1) +
  geom_smooth(aes(color = after_scale(fill)),
              method = "lm") +
  scale_fill_brewer(palette = "Dark2") +
  scale_alpha_continuous(trans = "log10",
                         range = c(0.05,.75),
                         n.breaks = 8) +
  coord_cartesian(expand = F) +
  facet_wrap(~Index, scales = "free") +
  baseTheme +
  theme(aspect.ratio = 1)

tibble(models = finalREModels %>% lapply(function(x) x$mer)) %>% 
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
    ),
    identity = map(data, ~ ggplot(.x, aes(fitted, fitted + residual, color = type)) +
                     geom_point(size = .1) +
                     geom_smooth(method = "lm") +
                     scale_color_brewer(palette = "Dark2") +
                     baseTheme +
                     theme(aspect.ratio = 1) +
                     facet_wrap(~Index, scales = "free",ncol=1)),
    predictors = map(data, ~ .x %>% 
                       pivot_longer(c(Richness,ShannonEvenness),names_to="Predictor") %>% 
                       ggplot(aes(value,residual,fill=type,color=type)) +
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
                       facet_grid(rows = vars(Index), cols = vars(Predictor), scales = "free", switch = "x") +
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
                       labs(x = NULL, y = "Residual",
                            alpha = "Density", fill = "Analysis\nMode")
    )),
    names = names(plots),
    width = c(4,3.5,4,4,4,4),
    height = c(2,3,2.5,2,4,4)
  ) %>% 
  mutate(save = pmap(list(names,plots,width,height), function(type,plot,w,h) {
    ggsave(paste0(type,"_REResidualPlot.png"),plot,
           type = "cairo-png", dpi = 400, scale = 3,
           width = w, height = h)
  }))

labelPT <- function(p,t,precision=2) map2(p,t, function(x,y) {
  xy <- sapply(c(x,y), function(z) {
    if (z == 0) {
      out <- "0"
    }
    else if (abs(z) > 0.01) {
      out <- str_extract(z,paste0("^.+\\..{0,",precision,"}"))
    } 
    else {
      out <- scales::scientific(z,precision) %>% 
        str_replace("e"," \\\\times \\\\, 10^{")
      if (str_detect(out,"[{]")) out <- paste0(out,"}")
    }
    out
  })
  
  paste0("$\\overset{",
         "P-value = ", xy[1],"^{",c(paste0(rep("*",3),collapse=""),paste0(rep("*",2),collapse=""),"*",".","")[which(x < c(0.0005, 0.005, 0.05, 0.1, 1))[1]],"}",
         "}{",
         "t-value = ", xy[2],
         "}$",
         collapse = "") %>% 
    latex2exp::TeX(output = "character")
})

# labelPT(0.0001,2.3)

finalREResultPlotData <- finalREResult %>% 
  filter(!str_detect(term,"Intercept")) %>% 
  # mutate(term = str_replace(term, "ShannonEvenness","Shannon Evenness") %>% 
  #          str_replace(":"," : "))
  mutate(term = factor(str_replace(term, "ShannonEvenness","Shannon Evenness") %>% 
                         str_replace(":"," : "),
                       levels = c(
                         "Bias",
                         "Richness",
                         "Shannon Evenness",
                         "Bias : Richness",
                         "Bias : Shannon Evenness"
                       ) %>% rev)) %>% 
  rename(estimate = 2, std.error = 3, t.value = 4, p.value = 5) %>% 
  group_by(Index) %>% 
  mutate(
    offset = max(abs(diff(range(c(estimate+std.error,estimate-std.error))))),
    direction = ifelse(estimate+std.error > max(estimate + std.error)*0.65, -1, 1)
  ) %>% 
  ungroup %>% 
  mutate(label = labelPT(p.value,t.value))

finalREResultPlot <- finalREResultPlotData %>% 
  ggplot(aes(estimate,Index,xmin=estimate-std.error*1.96,xmax=estimate+std.error*1.96,color=term)) +
  geom_vline(xintercept = 0, color = "black", size = .75) +
  geom_errorbar(width = .5,
                size = 1,
                position = position_dodge(1),
                show.legend = F) +
  geom_point(size = 4,
             position = position_dodge(1)) +
  geom_label(aes(x = estimate + direction * (1.96*std.error + 0.025*offset),
                 label = label),
             data = finalREResultPlotData %>% 
               mutate(estimate = ifelse(direction == -1,NA,estimate)),
             position = position_dodge(1),
             hjust = 0,
             parse = T,
             show.legend = F) +
  geom_label(aes(x = estimate + direction * (1.96*std.error + 0.025*offset),
                 label = label),
             data = finalREResultPlotData %>% 
               mutate(estimate = ifelse(direction == 1,NA,estimate)),
             position = position_dodge(1),
             hjust = 1,
             parse = T,
             show.legend = F) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(expand = expansion(.15)) +
  scale_y_discrete(expand = expansion()) +
  # coord_cartesian(ylim = c(.95,1.25)) +
  guides(color = guide_legend(override.aes = list(shape = 16,
                                                  size = 10,
                                                  alpha = 1),
                              ncol = 1,
                              reverse = T)) +
  labs(y = NULL, x = latex2exp::TeX("$\\hat{\\beta}$"), color = NULL) +
  ggpubr::theme_pubr() +
  theme(text = element_text(family = "CMU Serif"),
        title = element_text(face = "bold",
                             size = 20),
        legend.text = element_text(size = 20),
        strip.text = element_text(face = "bold",
                                  size = 20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = .4,
        legend.position = "right") +
  facet_wrap(~Index, scales = "free", ncol = 1, strip.position = "bottom")

ggsave("finalREResult.png",finalREResultPlot,
       type = "cairo-png", dpi = 400, height = 4, width = 4, scale = 3)


finalModelRE_predict <- finalModelData %>%
  rename("ShannonEvenness" = `Shannon\nEvenness`) %>%
  select(!c(FDis:FRic,fid,L2_KEY)) %>%
  summarize(across(where(is.numeric),~runif(20000, min(.x),max(.x))),
            across(c(Ecoregion),~sample(unique(.x),20000,T))) %>%
  mutate(rid = row_number()) %>%
  rowwise() %>%
  summarize(
    across(everything(),~rep(.x,2)),
    type = factor(c("Abundance","Presence/Absence"))
  ) %>%
  mutate(FDis = mgcv::predict.bam(finalREModels$FDis,.),
         FDiv = mgcv::predict.bam(finalREModels$FDiv,.),
         FEve = mgcv::predict.bam(finalREModels$FEve,.)) %>%
  group_by(rid) %>%
  summarize(
    across(FDis:FEve, ~(.x[1] - .x[2])/.x[1]),
    across(c(Richness,ShannonEvenness),first)
  )

finalModelRE_predict %>%
  pivot_longer(FDis:FEve,names_to="Index") %>%
  ggplot(aes(Richness,ShannonEvenness,z=value)) +
  stat_summary_2d(aes(color = after_scale(fill))) +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F) +
  facet_wrap(~Index)