library(tidyverse)
library(extrafont)

dat <- tibble(mu = sample(seq(0,50,5),6,T),
              sd = runif(6,5,20))
dat

nObs <- 15

sampleDat <- dat %>% 
  group_by(grp = row_number() %>% factor) %>% 
  summarize(sample = rnorm(nObs,mu,sd),
            ground = mu) %>% 
  ungroup

sampleDat %>% 
  ggplot(aes(grp,sample,group=grp)) +
  geom_violin(scale = "width") +
  stat_summary(fun.data = mean_cl_normal,
               color = "firebrick") +
  geom_point(position = "jitter") +
  geom_point(aes(y = ground),
             fill = "royalblue",
             shape = 21,
             size = 4) 

library(nlme)
library(multcomp)

library(rstatix)

model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)
}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)
}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)
}
 
Reduce(bind_rows,list(gls(sample ~ grp, data = sampleDat, correlation = varIdent(form = ~ 1 | grp)) %>% 
  glht(linfct = mcp(grp="Tukey")) %>% 
  summary %>% 
  broom::tidy() %>% 
  select(contrast,adj.p.value) %>% 
  rename(p.adj = 2) %>% 
  mutate(type = "GLS Tukey"),
TukeyHSD(aov(sample ~ grp, data = sampleDat)) %>% 
  broom::tidy() %>% 
  select(contrast,adj.p.value) %>% 
  mutate(contrast = str_replace(contrast, "-", " - ")) %>% 
  rename(p.adj = 2) %>% 
  mutate(type = "Regular Tukey"),
dunn_test(sampleDat, formula = sample ~ grp) %>% 
  mutate(contrast = paste0(group1, " - ", group2)) %>% 
  select(contrast,p.adj) %>% 
  mutate(type = "Dunn"),
pairwise.t.test(sampleDat$sample,sampleDat$grp,pool.sd = F, p.adjust.method = "bonferroni") %>% 
  broom::tidy() %>% 
  mutate(contrast = paste0(group1," - ",group2)) %>% 
  select(contrast,p.value) %>%
  rename(p.adj = 2) %>% 
  mutate(type = "T-test w.o.\npooled sd"),
apply(dat, 1, function(x) apply(dat, 1, function(y) {
  sapply(1:1000, function(DUMMY) mean(rnorm(nObs,x[1],x[2])) -  mean(rnorm(nObs,y[1],y[2]))) %>% 
    {
      min(1,(1 - abs(sum(sign(.)))/length(.))) 
    }
})) %>% 
  as_tibble %>% 
  mutate(group1 = row_number() %>% as.character()) %>% 
  pivot_longer(!group1,names_to="group2",values_to="p.value",names_prefix="V") %>%
  filter(group1 < as.numeric(group2)) %>% 
  mutate(contrast = paste0(group1," - ",group2)) %>% 
  select(contrast,p.value) %>% 
  rename(p.adj = 2) %>% 
  mutate(type = "Baseline"))) %>% 
  mutate(contrast = map_chr(contrast, function(x) str_split(x," - ") %>% 
                              unlist %>% 
                              as.numeric %>% 
                              sort %>% 
                              paste0(collapse = " - "))) %>% 
  # group_by(contrast,type) %>% 
  # summarize(
  #   across(everything(),~rep(.x,2)),
  #   contrast = c(contrast,
  #                str_split(contrast, " - ") %>% 
  #                  unlist %>% 
  #                  rev %>% 
  #                  paste0(collapse = " - "))
  # ) %>% 
  # ungroup %>% 
  # filter(type %in% c("Baseline","GLS Tukey", "Regular Tukey")) %>%
  ggplot(aes(contrast,p.adj,fill=type,group=type)) +
  geom_col(size = .5,
           color = "black",
           width = .5,
           position = position_dodge(.5)) +
  geom_hline(yintercept = 0.05) +
  scale_fill_brewer(palette = "Dark2") +
  baseTheme +
  coord_cartesian(expand = F) 
  # facet_wrap(~str_extract(contrast, "."),
  #            scales = "free_x")



bonfPLT <- pairwise.t.test(sampleDat$sample,sampleDat$grp,pool.sd = F, p.adjust.method = "bonferroni")$p.value %>% 
  t %>%  
  as.data.frame() %>% 
  rownames_to_column("group1") %>% 
  as_tibble %>%
  pivot_longer(!group1,names_to="group2",values_to="p.value") %>% 
  drop_na() %>% 
  mutate(across(1:2,~factor(.x,levels=unique(unlist(across(1:2)))))) %>%
  ggplot(aes(group2,group1,fill=p.value)) +
  geom_raster() +
  scale_fill_gradient2(na.value = "red",
                       low = "#FF3333",
                       high = "blue",
                       midpoint = log10(.05),
                       trans = "log10",
                       limits = c(10^-4,1),
                       labels = scales::label_pvalue()) +
  baseTheme +
  coord_equal(expand = F) +
  labs(x = NULL, y = NULL, title = "Pairwise t-test\nw. bonferroni")

holmPLT <- pairwise.t.test(sampleDat$sample,sampleDat$grp,pool.sd = F, p.adjust.method = "holm")$p.value %>% 
  t %>%  
  as.data.frame() %>% 
  rownames_to_column("group1") %>% 
  as_tibble %>%
  pivot_longer(!group1,names_to="group2",values_to="p.value") %>% 
  drop_na() %>% 
  mutate(across(1:2,~factor(.x,levels=unique(unlist(across(1:2)))))) %>%
  ggplot(aes(group2,group1,fill=p.value)) +
  geom_raster() +
  scale_fill_gradient2(na.value = "red",
                       low = "#FF3333",
                       high = "blue",
                       midpoint = log10(.05),
                       trans = "log10",
                       limits = c(10^-4,1),
                       labels = scales::label_pvalue()) +
  baseTheme +
  coord_equal(expand = F) +
  labs(x = NULL, y = NULL, title = "Pairwise t-test\nw. holm")


bhPLT <- pairwise.t.test(sampleDat$sample,sampleDat$grp,pool.sd = F, p.adjust.method = "BH")$p.value %>% 
  t %>%  
  as.data.frame() %>% 
  rownames_to_column("group1") %>% 
  as_tibble %>%
  pivot_longer(!group1,names_to="group2",values_to="p.value") %>% 
  drop_na() %>% 
  mutate(across(1:2,~factor(.x,levels=unique(unlist(across(1:2)))))) %>%
  ggplot(aes(group2,group1,fill=p.value)) +
  geom_raster() +
  scale_fill_gradient2(na.value = "red",
                       low = "#FF3333",
                       high = "blue",
                       midpoint = log10(.05),
                       trans = "log10",
                       limits = c(10^-4,1),
                       labels = scales::label_pvalue()) +
  baseTheme +
  coord_equal(expand = F) +
  labs(x = NULL, y = NULL, title = "Pairwise t-test\nw. BH")


tukeyPLT <- TukeyHSD(aov(sampleDat$sample~factor(sampleDat$grp)))$factor %>%
  as.data.frame() %>% 
  rownames_to_column("contrast") %>% 
  as_tibble() %>% 
  mutate(group = map(contrast, function(x) str_split(x,"-") %>% 
                       unlist %>% 
                       as.numeric %>% 
                       sort %>% 
                       setNames(1:2))) %>% 
  unnest_wider(group,names_sep="") %>% 
  relocate(contains("group")) %>% 
  select(!contrast) %>% 
  mutate(across(1:2,~factor(.x,levels=unique(unlist(across(1:2)))))) %>%
  ggplot(aes(group2,group1,fill=`p adj`)) +
  geom_raster() +
  scale_fill_gradient2(na.value = "red",
                       low = "#FF3333",
                       high = "blue",
                       midpoint = log10(.05),
                       trans = "log10",
                       limits = c(10^-4,1),
                       labels = scales::label_pvalue()) +
  baseTheme +
  coord_equal(expand = F) +
  labs(x = NULL, y = NULL, title = "Tukey Honest Signif-\nicance Differences",
       fill = "p.value")

simPLT <- apply(dat, 1, function(x) apply(dat, 1, function(y) {
  sapply(1:1000, function(DUMMY) mean(rnorm(nObs,x[1],x[2])^2) -  mean(rnorm(nObs,y[1],y[2])^2)) %>% 
    {
      1 - abs(sum(sign(.)))/length(.)
    }
})) %>% 
  as_tibble %>% 
  mutate(group1 = row_number() %>% as.character()) %>% 
  pivot_longer(!group1,names_to="group2",values_to="p.value",names_prefix="V") %>% 
  mutate(across(1:2,~factor(.x,levels=1:10))) %>%
  filter(as.numeric(as.character(group1)) > as.numeric(as.character(group2))) %>%
  ggplot(aes(group1,group2,fill=p.value)) +
  geom_raster() +
  scale_fill_gradient2(na.value = "red",
                       low = "#FF3333",
                       high = "blue",
                       midpoint = log10(.05),
                       trans = "log10",
                       limits = c(10^-4,1),
                       labels = scales::label_pvalue()) +
  baseTheme +
  coord_equal(expand = F) +
  labs(x = NULL, y = NULL, title = "Simulation")

library(patchwork)
wrap_plots(list(bonfPLT,holmPLT,bhPLT,simPLT,tukeyPLT),
           guides = "collect")


