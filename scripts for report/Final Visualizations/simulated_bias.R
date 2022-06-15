library(asgerbachelor)
library(tidyverse)
library(magrittr)

minSP <- 10
maxSP <- 100
nTr <- 50

n <- 10000000

communities <- rlnorm(n,-4,2) %>%
  raise_to_power(-1) %>% 
  matrix(ncol=maxSP) %>%
  apply(1, function(x) replace(x,sample(maxSP,maxSP-sample(minSP:maxSP,1),F),0)) %>% 
  t %>% 
  divide_by(rowSums(.)) 

traits <- matrix(rnorm(maxSP*nTr),ncol=nTr)

int_to_vec <- function(x) {
  x %>% 
    str_extract_all("[[:digit:][\\.]-]+") %>% 
    unlist %>% 
    as.numeric %>% 
    set_names(c("min","max"))
}

rwDat <- tibble(
  Richness = rowSums(communities>0),
  Evenness = apply(communities,1,function(x) -sum(x*log(x),na.rm=T))/log(Richness),
  Bias = functional_dispersion(traits,a=communities) - functional_dispersion(traits,a=communities>0)
) %>% 
  dplyr::filter(Evenness != 0) %>% 
  group_by(Evenness = cut_interval(Evenness, 30, boundary = 0, dig.lab = 10),
           Richness = cut_interval(Richness, 30, boundary = 0, dig.lab = 10)) %>% 
  summarize(Bias = abs(mean(Bias)),
            .groups = "drop") %>% 
  mutate(across(where(is.factor), ~map(.x,int_to_vec))) %>% 
  unnest_wider(where(is.list),names_sep="_")

randomWeightFDisPlot <- rwDat %>% 
  ggplot(aes(fill=Bias)) +
  geom_rect(aes(xmin = Evenness_min, xmax = Evenness_max,
                ymin = Richness_min, ymax = Richness_max)) + 
  scale_fill_distiller(palette = "Reds",
                       direction = 1,
                       limits = range(rwDat$Bias),
                       breaks = range(rwDat$Bias),
                       labels = c("Low","High")) +
  # scale_fill_distiller(palette = "Spectral",
  #                      limits = c(-1,1)*max(abs(range(rwDat$diff)))) +
  scale_x_continuous(breaks = seq(0.05,1,.2)) +
  scale_y_continuous(breaks = seq(10,100,30)) +
  coord_cartesian(expand = F) +
  ggpubr::theme_pubr() +
  labs(x = "Shannon Evenness",
       y = "Species Richness",
       fill = "Bias") +
  baseTheme +
  theme(legend.position = "right",
        legend.key.height = unit(3,"lines"),
        legend.text = element_text(size = 14),
        title = element_text(size = 20),
        plot.title = element_text(hjust = .5,
                                  size = 24),
        aspect.ratio = 1)

ggsave("exampleRandomWeights.png",randomWeightFDisPlot,
       type = "cairo-png",width = 3, height = 2.4, scale = 2)
