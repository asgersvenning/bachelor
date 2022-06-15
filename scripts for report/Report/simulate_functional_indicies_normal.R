library(magrittr)
library(tidyverse)
library(asgerbachelor)

s <- 200
n <- 500

virtual_normal_species <- rnorm(s*2) %>% 
  matrix(ncol = 2)

virtual_lognormal_communities <- rlnorm(s*n,
                                        -4,
                                        # rep(rnorm(n,-2,1),each=s),
                                        1.25
                                        # rep(runif(n,.75,1.5),each=s)
) %>%
  raise_to_power(-1) %>%
  matrix(ncol=s) %>%
  apply(1, function(x) replace(x,sample(s,s-sample(5:(s-50),1),F),0)) %>% 
  t %>% 
  divide_by(rowSums(.)) %>% 
  set_colnames(paste0("SP",1:s))

tibble(rich = rowSums(virtual_lognormal_communities>0),
       even = apply(virtual_lognormal_communities,1,function(x) -sum(x*log(x),na.rm=T))/log(rich)) %>% 
  ggplot(aes(rich,even)) +
  geom_bin_2d(aes(color = after_scale(fill)),
              bins = 15) +
  scale_fill_viridis_c(trans = "log10") +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(log(x)),
              method.args = list(famiy = "quasibinomial"))


virtual_communities_pa_and_ab <- rbind(virtual_lognormal_communities,virtual_lognormal_communities>0)

functional_indices_normal <- tibble(
  richness = rowSums(virtual_lognormal_communities>0),
  evenness = apply(virtual_lognormal_communities, 1, 
                   function(x) -sum(x/sum(x) * log(x/sum(x)), na.rm = T))/log(richness),
  FDis = functional_dispersion(virtual_normal_species,a=virtual_lognormal_communities) - 
    functional_dispersion(virtual_normal_species,a=virtual_lognormal_communities>0),
  FEve = functional_evenness(virtual_normal_species,a=virtual_lognormal_communities) - 
    functional_evenness(virtual_normal_species,a=virtual_lognormal_communities>0),
  FDiv = functional_divergence(virtual_normal_species,a=virtual_lognormal_communities) -
    functional_divergence(virtual_normal_species,a=virtual_lognormal_communities>0)
) 

# functional_indices_normal_plot <- functional_indices_normal %>% 
#   pivot_longer(!c(richness,evenness, type)) %>% 
#   ggplot(aes(richness,value,color=name)) +
#   geom_point() +
#   stat_summary_bin(color = "black", bins = 10) +
#   scale_color_brewer(palette = "Dark2") + 
#   facet_grid(rows = vars(name), cols = vars(type), scales = "free_y") +
#   guides(color = guide_legend(override.aes = list(shape = 15,
#                                                   size = 10))) + 
#   ggpubr::theme_pubr() +
#   theme(text = element_text(family = "CMU Serif"),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(hjust = .5),
#         strip.text = element_text(face = "bold",
#                                   size = 18),
#         title = element_text(face = "bold",
#                              size = 18),
#         legend.position = "right",
#         plot.title = element_text(hjust = .5, size = 20, face = "plain")) +
#   labs(x = "Species Richness", y = "Index Value", color = "Functional\nIndex",
#        title = "Normal Virtual Species")
# 
# ggsave("functional_indices_normal.png", functional_indices_normal_plot,
#        type = "cairo-png", dpi = 300, width = 5, height = 6, scale = 2)

library(lmerTest)

tfin <- functional_indices_normal %>% 
  mutate(logmean = virtual_lognormal_communities %>% 
           log %>% 
           apply(1, function(x) mean(x[is.finite(x)])),
         shannon = virtual_lognormal_communities %>% 
           apply(1, function(x) sum(x*log(x),na.rm=T)))

tfin %>% 
  mutate(richness = log(richness)) %>% 
  pivot_longer(c(FDis,FDiv,FEve),names_to="Index") %>% 
  pivot_longer(c(richness,evenness),names_to="Structure",values_to="predictor") %>% 
  ggplot(aes(predictor,value,color=Index)) +
  geom_point(shape=21,fill="gray75") +
  geom_smooth(method="lm", formula = y ~ x) +
  coord_cartesian(expand = F) +
  facet_wrap(~Structure,scales="free") +
  theme_bw() +
  theme(aspect.ratio = 1) 

gam(FEve ~ s(log(richness)) + s(evenness),data=tfin) %>% 
  broom::tidy()

lm(FEve ~ evenness + richness, data = tfin) %>% 
  anova %>% 
  as.data.frame %>% 
  rownames_to_column("Effect") %>% 
  as_tibble

functional_indices_normal_lm <- functional_indices_normal %>% 
  mutate(richness = log(richness),
         # evenness = log(evenness/(1-evenness)),
         across(c(richness,evenness),~as.vector(scale(.x))),
         shannon = virtual_lognormal_communities %>% 
           apply(1, function(x) sum(x*log(x),na.rm=T))) %>% 
  pivot_longer(c(FDis,FDiv,FEve)) %>% 
  nest(dat = !name) %>% 
  mutate(
    Without = map(dat, ~car::Anova(lm(value ~ richness, data = .x),type="2") %>%
                as.data.frame %>% 
                rownames_to_column("Effect") %>% 
                as_tibble),
    With = map(dat, ~car::Anova(lm(value ~ evenness + richness, data = .x),type="2") %>%
                as.data.frame %>% 
                rownames_to_column("Effect") %>% 
                as_tibble)
  ) %>% 
  select(name,With,Without) %>% 
  pivot_longer(c(With,Without),names_to="type") %>% 
  unnest(where(is.list))

functional_indices_normal_lm %>% 
  mutate(Effect = ifelse(Effect == "Residuals","Residual",Effect) %>% 
           factor(levels = c("richness","evenness","Residual"))) %>%
  arrange(name,Effect,type) %>% 
  select(1:4,6:7) %>% 
  filter(Effect == "richness") %>% 
  select(!Effect)

write_csv(functional_indices_normal_lm,"functional_indices_normal_lm.csv")
