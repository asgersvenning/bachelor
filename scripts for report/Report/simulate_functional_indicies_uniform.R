library(magrittr)
library(tidyverse)
library(asgerbachelor)

virtual_uniform_species <- runif(100*5) %>% 
  matrix(ncol = 5)

virtual_lognormal_communities <- sapply(rep(10:100,3), function(x) {
  s <- sample(100,x)
  w <- rlnorm(x,-4,1.25) %>%
    divide_by(sum(.))
  w <- rep(0,100) %>% 
    inset(s,w)
  
  set_names(w,paste0("SP",1:length(w)))
}) %>% 
  t

virtual_communities_pa_and_ab <- rbind(virtual_lognormal_communities,virtual_lognormal_communities>0)

functional_indices_uniform <- tibble(
  richness = rep(rep(10:100,3),2),
  FDis = functional_dispersion(virtual_uniform_species,,virtual_communities_pa_and_ab),
  FEve = functional_evenness(virtual_uniform_species,,virtual_communities_pa_and_ab),
  FRic = functional_richness(virtual_uniform_species,,virtual_communities_pa_and_ab),
  FDiv = functional_divergence(virtual_uniform_species,,virtual_communities_pa_and_ab)
) %>% 
  mutate(type = rep(c("Abundance","Presence/Absence"),each=n()/2))

functional_indices_uniform_plot <- functional_indices_uniform %>% 
  pivot_longer(!c(richness,type)) %>% 
  ggplot(aes(richness,value,color=name)) +
  geom_point() +
  stat_summary_bin(color = "black", bins = 10) +
  scale_color_brewer(palette = "Dark2") + 
  facet_grid(rows = vars(name), cols = vars(type), scales = "free_y") +
  guides(color = guide_legend(override.aes = list(shape = 15,
                                                  size = 10))) + 
  ggpubr::theme_pubr() +
  theme(text = element_text(family = "CMU Serif"),
        legend.text = element_text(size = 12),
        legend.title = element_text(hjust = .5),
        strip.text = element_text(face = "bold",
                                  size = 18),
        title = element_text(face = "bold",
                             size = 18),
        legend.position = "right",
        plot.title = element_text(hjust = .5, size = 20, face = "plain")) +
  labs(x = "Species Richness", y = "Index Value", color = "Functional\nIndex",
       title = "Uniform Virtual Species")

ggsave("functional_indices_uniform.png", functional_indices_uniform_plot,
       type = "cairo-png", dpi = 300, width = 5, height = 6, scale = 2)

library(kableExtra)
functional_indices_uniform_lm <- functional_indices_uniform %>% 
  pivot_longer(!c(richness,type)) %>% 
  nest(dat = !name) %>% 
  mutate(mod = map(dat, ~lm(value ~ richness*type, data = .x)),
         res = map(mod, ~.x %>% 
                     summary %>% 
                     coefficients %>% 
                     as.data.frame %>% 
                     rownames_to_column("Effect") %>% 
                     as_tibble)) %>% 
  select(name,res) %>% 
  unnest(res) 

write_csv(functional_indices_uniform_lm,"functional_indices_uniform_lm.csv")
