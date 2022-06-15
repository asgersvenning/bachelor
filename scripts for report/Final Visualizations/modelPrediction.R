library(patchwork)

predictedIndices <- finalModelData %>% 
  filter(type == "Abundance") %>% 
  # mutate(nPlots = 100) %>% 
  select(type,Richness,ShannonEvenness,nPlots,Ecoregion,x,y) %>% 
  summarize(across(c(Richness,ShannonEvenness),~runif(10000,min(.x),max(.x))),
            wID = sample(n(),10000,T),
            x = x[wID],
            y = y[wID],
            nPlots = sample(nPlots,10000,T),
            type = sample(c("Abundance","Presence/Absence"),10000,T),
            .groups = "drop") %>% 
  mutate(FEve = predict(finalModels %>% 
                          filter(form == "type*(Richness + ShannonEvenness)" &
                                   Index == "FEve") %>% 
                          pull(model) %>% {.[[1]]$gam},.),
         FDis = predict(finalModels %>% 
                          filter(form == "type*(Richness + ShannonEvenness)" &
                                   Index == "FDis") %>% 
                          pull(model) %>% {.[[1]]$gam},.),
         FDiv = predict(finalModels %>% 
                          filter(form == "type*(Richness + ShannonEvenness)" &
                                   Index == "FDiv") %>% 
                          pull(model) %>% {.[[1]]$gam},.)) %>%
  mutate(
    Richness = exp((0.837*Richness+3.17906)),
    ShannonEvenness = ShannonEvenness*0.0937 + 0.8026
  ) 

predictedIndices %>% 
  mutate(Richness = log(Richness)) %>% 
  pivot_longer(c(FEve,FDis,FDiv),names_to="Index",values_to="IV") %>%
  pivot_longer(c(Richness,ShannonEvenness)) %>% 
  nest(dat = !type) %>% 
  mutate(plt = map2(dat,type, function(a,b) {
    ggplot(a,aes(value,IV)) +
      geom_hex(aes(color=after_scale(fill))) +
      stat_summary_bin(color="firebrick1") +
      geom_smooth() +
      scale_fill_viridis_c(trans = "log10") +
      facet_grid(rows = vars(Index), cols = vars(name), scales = "free") +
      coord_cartesian(expand = F) +
      baseTheme +
      theme(aspect.ratio = 1,
            plot.title = element_text(hjust = .5)) +
      labs(title = b, x = NULL, y = NULL)
  })) %>% 
  pull(plt) %>% 
  wrap_plots()


predictedIndices %>% 
  pivot_longer(c(FDis,FDiv,FEve)) %>% 
  nest(dat = !name) %>% 
  mutate(plt = map2(dat,name, function(d,n) {
    d %>% 
      ggplot(aes(Richness,ShannonEvenness,z=value)) +
      stat_summary_2d(aes(color=after_scale(fill)),
                       bins = 15,
                      drop = F) +
      scale_fill_viridis_c() +
      scale_x_log10() +
      coord_cartesian(expand = F) +
      baseTheme +
      theme(aspect.ratio = 1) +
      labs(fill = n)
  })) %>% 
  pull(plt) %>% 
  patchwork::wrap_plots()

modelDataRaw$Richness %>% log %>% scale
modelDataRaw$`Shannon\nEvenness` %>% scale