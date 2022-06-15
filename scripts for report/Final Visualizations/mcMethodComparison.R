library(tidyverse)
library(extrafont)
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

model.matrix.gam <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)
}
model.frame.gam <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)
}
terms.gam <- function(object, ...) {
  terms(model.frame(object), ...)
}

res <- Reduce(bind_rows,lapply(1:10, function(dummy) {
  dat <- tibble(mu = sample(0:10 * 5,5,T),
                sd = runif(5,1,5))
  
  nObs <- 25
  
  sampleDat <- dat %>% 
    group_by(grp = row_number() %>% factor) %>% 
    summarize(sample = rnorm(nObs,mu,sd),
              ground = mu) %>% 
    ungroup %>% 
    mutate(grp = factor(grp))
  
  Reduce(bind_rows,list(gls(sample ~ grp, data = sampleDat, 
                            correlation = varIdent(form = ~ 1 | grp)) %>% 
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
                        dunn_test(sampleDat, formula = sample ~ grp, p.adjust.method = "BH") %>% 
                          mutate(contrast = paste0(group1, " - ", group2)) %>% 
                          select(contrast,p.adj) %>% 
                          mutate(type = "Dunn (BH)"),
                        dunn_test(sampleDat, formula = sample ~ grp, p.adjust.method = "holm") %>% 
                          mutate(contrast = paste0(group1, " - ", group2)) %>% 
                          select(contrast,p.adj) %>% 
                          mutate(type = "Dunn (holm)"),
                        pairwise.t.test(sampleDat$sample,sampleDat$grp,pool.sd = F, p.adjust.method = "holm") %>% 
                          broom::tidy() %>% 
                          mutate(contrast = paste0(group1," - ",group2)) %>% 
                          select(contrast,p.value) %>%
                          rename(p.adj = 2) %>% 
                          mutate(type = "T-test w.o.\npooled sd (holm)"),
                        pairwise.t.test(sampleDat$sample,sampleDat$grp,pool.sd = F, p.adjust.method = "BH") %>% 
                          broom::tidy() %>% 
                          mutate(contrast = paste0(group1," - ",group2)) %>% 
                          select(contrast,p.value) %>%
                          rename(p.adj = 2) %>% 
                          mutate(type = "T-test w.o.\npooled sd (BH)"),
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
    group_by(contrast) %>% 
    mutate(type1 = ifelse(p.adj < p.adj[type == "Baseline"] & p.adj < 0.05, 1, 0),
           type2 = ifelse(p.adj > p.adj[type == "Baseline"] & p.adj > 0.05, 1, 0)) %>% 
    ungroup %>% 
    group_by(type) %>% 
    summarize(across(c(type1,type2),sum),
              n = n())}))

res %>% 
  group_by(type) %>% 
  summarize(across(everything(),sum)) %>% 
  mutate(none = n - type1 - type2) %>% 
  mutate(across(c(type1,type2,none),~.x/n)) %>% 
  pivot_longer(c(type1,type2,none),names_to="error") %>% 
  mutate(error = map_chr(error, ~switch(.x,
                                        "none" = "None",
                                        "type1" = "Type 1",
                                        "type2" = "Type 2"))) %>% 
  filter(type != "Baseline") %>% 
  ggplot(aes(type,value,fill=error)) +
  geom_col(color = "black", size = .75, width = .75,
           position = "dodge", key_glyph = draw_key_point) +
  coord_cartesian(expand = F) +
  baseTheme +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Statistical test", y = "Error Frequency", fill = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 10,
                                                 stroke = 1,
                                                 shape = 21))) +
  theme(legend.text = element_text(size = 14))
