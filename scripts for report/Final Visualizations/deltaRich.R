library(asgerbachelor)
library(tidyverse)
library(extrafont)
library(mgcv)
library(lme4)
library(Hmisc)

dRich <- FIA %>%
  filter(ID %in% tSP$ID) %>%
  group_by(ID) %>%
  mutate(
    samples = sapply(unique(INVYR), function(yr) {
      mean(samples[yr==INVYR],na.rm=T)
    }) %>%
      sum
  ) %>%
  group_by(ID,INVYR) %>% 
  summarise(
    samples = first(samples),
    Richness = n(),
    .groups = "drop"
  ) %>% 
  # mutate(eco = communityCoords[sapply(ID, function(x) which(x == sort(unique(ID)))),]) %>% 
  # unnest(eco) %>% 
  group_by(ID) %>% 
  arrange(desc(INVYR)) %>% 
  summarise(
    deltaRich = weighted.mean((Richness/samples-lag(Richness/samples))/(INVYR-lag(INVYR)),w=(samples + lag(samples))/2,na.rm=T),
    baselineRich = first(Richness/samples),
    samples = sum(samples)
  )
  # mutate(INVYR = INVYR - min(INVYR)) %>% 
  # summarize(
  #   estimates = glm(Richness ~ INVYR, weights = samples/sum(samples), family = "quasipoisson") %>% 
  #     coefficients(),
  #   terms = names(estimates)
  # )

# gamm_Mod <- mgcv::bam(Richness ~ INVYR*L2_KEY + s(x, y) + s(samples),
#                       data = dRich,
#                       weights = samples,
#                       family = "quasipoisson")
# 
# dRich %>% mutate(fitted = gamm_Mod %>% fitted) %>% 
#   ggplot(aes(Richness,fitted)) +
#   stat_bin_2d() +
#   geom_smooth() + 
#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") + 
#   scale_fill_viridis_c(trans = "log10") +
#   scale_x_log10() +
#   scale_y_log10()
# 
# dRich %>% mutate(fitted = gamm_Mod %>% fitted) %>% 
#   ggplot(aes(x,y,fill=sqrt(abs(log(Richness/fitted))))) +
#   geom_raster() +
#   scale_fill_viridis_c() +
#   coord_equal()
# 
# par(mfrow = c(2,2))
# gamm_Mod %>% 
#   mgcv::gam.check()
# par(mfrow = c(1,1))
# 
# gamm_Mod %>% 
#   summary()

allIndices %>% 
  rowwise %>% 
  mutate(dRich = dRich %>% 
           # pivot_wider(ID,names_from=terms,values_from=estimates) %>% 
           # mutate(deltaRich = exp(INVYR)/exp(`(Intercept)`)) %>% 
           bind_cols(communityCoords) %>% 
           list) %>% 
  unnest(everything()) %>% 
  pivot_longer(FDis:FRic,names_to="Index") %>% 
  filter(Index != "FRic") %>% 
  ggplot(aes(value,deltaRich,fill=type)) +
  geom_point(shape = 21, alpha = .25) +
  geom_smooth() +
  facet_grid(cols = vars(Index), rows = vars(type), scales = "free") +
  scale_fill_brewer(palette = "Dark2") + 
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major = element_line(color = "gray65"),
        text = element_text(family = "CMU Serif")) 



lmDat <- allIndices %>% 
  rowwise %>% 
  mutate(dRich = dRich %>% 
           # pivot_wider(ID,names_from=terms,values_from=estimates) %>% 
           # mutate(deltaRich = exp(INVYR)/exp(`(Intercept)`)) %>% 
           bind_cols(communityCoords) %>% 
           list,
         Richness = rowSums(fSP>0) %>% list,
         ShEve = (fSP %>% 
                    apply(1, 
                          function(x) -sum(x*log(x),na.rm=T)) / log(Richness)) 
         %>% list) %>% 
  unnest(everything()) 

library(lmerTest)
stepLM <- lmerTest::lmer(deltaRich ~ (FDis + FDiv + FEve + FRic)*(baselineRich + Richness + ShEve) + (1 | L2_KEY / fid),
             data = lmDat %>% 
               mutate(samples = samples/sum(samples)) %>% 
               # mutate(across(c(FDis,FDiv,FEve,FRic,baselineRich,Richness,ShEve,samples),~(.x-min(.x,na.rm=T))/diff(range(.x,.na.rm=T)))) %>% 
               group_by(type) %>% 
               mutate(fid = row_number() %>% factor),
               # ungroup %>% 
               # mutate(across(where(is.numeric), ~(.x-min(.x,na.rm=T))/diff(range(.x,na.rm=T)))),
             weight = samples
               # filter(type == "Abundance") %>%
               # dplyr::select(!type)
             ) %>% 
  lmerTest::step() %>% 
  get_model()

summary(stepLM) 
par(mfrow = c(2,2))
plot(stepLM)

stepLM %>% 
  broom::tidy() %>% 
  filter(!str_detect(term,"KEY") & p.value < .1) %>% 
  arrange(p.value)


lmDat %>% 
  # filter(type == "Abundance") %>% 
  ggplot(aes(baselineRich,deltaRich,z=FRic)) +
  stat_summary_hex(aes(color = after_stat(value)),
                   bins = 25) +
  # geom_point(size = 2,
             # alpha = .5) +
  scale_fill_viridis_c(
    # trans = "log10",
    option = "D"
  ) +
  scale_color_viridis_c(
    # trans = "log10",
    option = "D",
    guide = "none"
  ) +
  scale_shape_manual(values = c(21,23)) + 
  theme_dark() +
  facet_wrap(~type, scales = "free_y")



lmDat %>% 
  # filter(type == "Abundance") %>% 
  ggplot(aes(Richness,deltaRich,z=ShEve)) +
  stat_summary_hex(aes(color = after_stat(value)),
                   bins = 25) +
  # geom_point(size = 2,
  # alpha = .5) +
  scale_fill_viridis_c(
    # trans = "log10",
    option = "D"
  ) +
  scale_color_viridis_c(
    # trans = "log10",
    option = "D",
    guide = "none"
  ) +
  scale_shape_manual(values = c(21,23)) + 
  theme_dark() +
  facet_wrap(~type, scales = "free_y")
