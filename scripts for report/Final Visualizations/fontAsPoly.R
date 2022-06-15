library(gglogo)
library(tidyverse)
library(extrafont)

CMU_poly <- createPolygons(c(LETTERS,letters,0:9,".",",","!","?"),font="CMU Serif")

CMU_tidy <- CMU_poly %>% 
  as_tibble %>%
  select(x,y,group) %>% 
  nest(poly = !group) 

string_as_poly <- function(str,poly) {
  if (str_detect(str,paste0(c("[",poly$group,"[ ]]"),collapse=""),negate=T)) {
    warnings("Unknown characters detected")
    str <- str_extract_all(str,paste0("[",poly$group,"[ ][\n]]")) %>% 
      unlist %>% 
      paste0(collapse="")
  }
  str <- str_extract_all(str,"\n|.") %>% 
    unlist
  
  
  tibble(letter = str,
         group = 1:length(str)) %>% 
    mutate(line = ifelse(letter == "\n",1,0),
           line = cumsum(line)) %>% 
    filter(letter != "\n") %>% 
    mutate(nest = map(letter,function(x) {
      wL <- which(poly$group==x)
      
      if (length(wL) == 0) {
        tibble() 
      }
      else {
        poly$poly[[wL]] %>% 
          mutate(x = x - min(x))
      }
    })) %>% 
    mutate(width = map_dbl(nest, function(n) {
      if (all(dim(n) == c(0,0))) return(.25) else abs(diff(range(n$x)))
    })) %>% 
    group_by(line) %>% 
    mutate(x_off = cumsum(lag(width,default = 0)+0.1),
           y_off = -line) %>% 
    ungroup
}

string_as_poly(stringi::stri_rand_lipsum(1) %>% str_wrap(40),CMU_tidy) %>% 
  unnest(nest) %>% 
  ggplot(aes(x+x_off,y+y_off,group=group,fill=sqrt((x_off-mean(x_off))^2 + (y_off-mean(y_off))^2))) +
  geom_polygon(show.legend = F) + 
  scale_fill_viridis_c() +
  theme_void() +
  coord_equal() 


