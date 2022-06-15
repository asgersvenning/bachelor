n <- 1000000

tMat <- rlnorm(n,-4,1) %>%
  raise_to_power(-1) %>% 
  matrix(ncol=200) %>%
  apply(1, function(x) replace(x,sample(200,sample(0:190,1),F),0)) %>% 
  t %>% 
  # divide_by(rowSums(.)) %>%
  # apply(1, function(x) {
  #   rarTab <- sample(200,200,T,prob=x) %>% table
  #   vector(mode="integer",length=200) %>%
  #     replace(as.numeric(names(rarTab)),rarTab)
  # }) %>%
  # t %>%
  divide_by(rowSums(.))

tTr <- matrix(runif(200*5),ncol=5)

tDat <- tibble(richness = rowSums(tMat>0),
               evenness = apply(tMat,1,function(x) -sum(x*log(x),na.rm=T))/log(richness),
               bias = functional_dispersion(tTr,a=tMat) - functional_dispersion(tTr,a=tMat>0))

tDat %>% 
  ggplot(aes(richness,evenness, z=bias)) +
  stat_summary_2d(aes(color=after_scale(fill)),
                  bins = 15) +
  scale_fill_viridis_c() +
  stat_summary_bin(fun.data = mean_sdl,
                   fun.args = list(mult = .1)) +
  geom_smooth(method = "lm", formula = y ~ log(x)) 

tDat %>% 
  group_by(log_int = cut_interval(log(richness),50,dig.labs=10)) %>% 
  summarize(
    "Log-Variance" = var(evenness) %>% log,
    Mean = mean(evenness)
  ) %>% 
  mutate(richness = map(log_int, str_extract_all, pattern = "[[:digit:]\\.-]+") %>% 
           map(~.x %>% 
                 unlist %>% 
                 as.numeric %>% 
                 set_names(c("min","max")))
  ) %>% 
  unnest_wider(richness,names_sep="_") %>% 
  mutate(richness = rowMeans(across(contains("rich"))),
         across(contains("richness"),exp)) %>% 
  pivot_longer(2:3) %>% 
  ggplot(aes(richness,value,fill=name)) +
  geom_point(size = 3, stroke = 1,
             shape = 21) +
  geom_smooth(method = "lm",
              color = "black") +
  scale_x_log10() +
  facet_wrap(~name,scales="free") +
  coord_cartesian(expand = F) +
  baseTheme +
  theme(aspect.ratio = 1) +
  labs(x = "Species Richness", y = NULL)



modelData %>% 
  select(Richness,nPlots,`Shannon\nEvenness`) %>% 
  rename("Evenness" = 3) %>% 
  # mutate(Evenness = residuals(mgcv::gam(Evenness ~ s(log(nPlots),log(Richness))))) %>% 
  group_by(log_int = cut_number(log(Richness),25,dig.lab=20),
           nPlots = cut_number(nPlots,2,dig.lab=5))  %>% 
  summarize(
    "Variance" = var(Evenness),
    Mean = mean(Evenness),
    n = n(),
    .groups = "drop"
  ) %>% 
  drop_na %>% 
  mutate(Richness = map(log_int, str_extract_all, pattern = "[[:digit:]\\.-]+") %>% 
           map(~.x %>% 
                 unlist %>% 
                 as.numeric %>% 
                 set_names(c("min","max")))
  ) %>% 
  unnest_wider(Richness,names_sep="_") %>% 
  mutate(Richness = rowMeans(across(contains("Richness"))),
         across(contains("richness"),exp)) %>% 
  pivot_longer(3:4) %>% 
  ggplot(aes(Richness,value,fill=name,size=n,weight=n)) +
  geom_point(stroke = 1,
             shape = 21) +
  geom_smooth(color = "black",
              method="gam",
              show.legend = F) +
  scale_size_continuous(limits = c(0,150)) + 
  # scale_x_log10() +
  facet_grid(rows = vars(name), cols = vars(nPlots),scales="free_y") +
  coord_cartesian(expand = T) +
  baseTheme +
  theme(aspect.ratio = 1) +
  labs(x = "Species Richness", y = NULL)


modelData %>% 
  select(Richness,nPlots,`Shannon\nEvenness`) %>% 
  rename("Evenness" = 3) %>% 
  ggplot(aes(Richness,Evenness)) +
  geom_point() +
  stat_summary_bin(fun.data = mean_sdl,
                   color = "red",
                   bins = 10) + 
  geom_smooth() +
  # scale_x_log10() +
  facet_wrap(~cut_number(nPlots,8,dig.lab=5))

