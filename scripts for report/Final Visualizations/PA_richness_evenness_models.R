source("tidyData_and_calculate_indices.R")
library(mgcv)
library(gamm4)
library(ggtext)

baseTheme <-  ggpubr::theme_pubr() +
  theme(text = element_text(family = "CMU Serif"),
        title = element_text(face = "bold",
                             size = 14),
        strip.text = element_text(face = "bold",
                                  size = 14),
        legend.position = "right")

finalModelData <- modelDataRaw %>%  
  bind_cols(bind_rows(communityCoords,communityCoords)) %>% 
  mutate(Ecoregion = factor(Ecoregion),
         type = factor(type),
         Richness = log(Richness),
         across(c(Richness,`Shannon\nEvenness`), ~ scale(.x) %>% as.vector)) %>%
  rename("ShannonEvenness" = `Shannon\nEvenness`)

# finalModels <- tibble(form = c(
#   "1",
#   "type",
#   "type + Richness",
#   "type*Richness",
#   "type + ShannonEvenness",
#   "type*ShannonEvenness",
#   "type*(Richness + ShannonEvenness)"
# )) %>%
#   mutate(Index = list(c("FDis",
#                         "FDiv",
#                         "FEve"))) %>%
#   unnest(Index) %>%
#   mutate(model = map2(Index, form, function(x,y) {
#     mod <- paste0(x, " ~ ", y, " + s(x,y,log(nPlots), by = type)")
#     print(mod)
#     
#     gamm4(formula(mod),
#           random = ~ (1 | Ecoregion/fid),
#           data = finalModelData)
#   }))
# 
# write_rds(finalModels,"finalModels.rds")

finalModels <- read_rds("finalModels.rds")

finalModels %>% 
  mutate(tabName = map_chr(form, ~switch(.x,
                                         "1"="1",
                                         "type"="bias",
                                         "type + Richness"="bias_p_rich",
                                         "type*Richness"="bias_x_rich",
                                         "type + ShannonEvenness" = "bias_p_even",
                                         "type*ShannonEvenness"="bias_x_even",
                                         "type*(Richness + ShannonEvenness)"="bias_x_rich_even")),
         tab = map2(model, form, function(x,y) {
           out <- summary(x$gam)$p.table %>% 
             as.data.frame %>% 
             rownames_to_column %>% 
             as_tibble %>% 
             mutate(rowname = rowname %>%
                      str_remove("^X") %>% 
                      str_replace("typePresence/Absence","Bias")) %>% 
             rename(term = rowname) 
         }),
         AIC = map_dbl(model,~AIC(.x$mer))) %>% 
  select(Index,tabName,tab,AIC) %>% 
  unnest(tab) %>% 
  nest(tab = !c(Index,tabName,AIC)) %>% 
  write_rds("allTab.rds")


finalModelAIC_dat <- finalModels %>% 
  mutate(mer = map(model, ~.x$mer)) %>% 
  group_by(Index) %>% 
  summarize(
    an = with(mer %>% set_names(form),
              do.call(anova, 
                      c(lapply(form, as.name),refit=F))) %>% 
      as.data.frame %>% 
      rownames_to_column("model") %>% 
      as_tibble,
    r.sq = map_dbl(model,~summary(.x$gam)$r.sq)
  ) %>% 
  unnest(where(is.list)) %>% 
  # mutate(form = map_chr(model,function(x) {
  #   switch(x,
  #          "1" = "$\\beta$",
  #          "type" = "$\\beta + \\beta^*$",
  #          "type + Richness" = "$\\beta + \\beta^* + \\alpha_R$",
  #          "type*Richness" = "$\\beta + \\alpha_R + \\beta^* + \\alpha^*_R$",
  #          "type + ShannonEvenness" = "$\\beta + \\beta^* + \\alpha_{SE}$",
  #          "type*ShannonEvenness" = "$\\beta + \\alpha_{SE} + \\beta^* + \\alpha^*_{SE}$",
  #          "type*(Richness + ShannonEvenness)" = "$\\beta + \\alpha_R + \\alpha_{SE} + \\beta^* + \\alpha^*_R + \\alpha^*_{SE}$") %>% latex2exp::TeX(out = "character")
  # })) %>%
  mutate(form = map_chr(model,function(x) {
    switch(x,
           "1" = "$\\beta$",
           "type" = "$\\beta + \\beta^{Bias}$",
           "type + Richness" = "$\\beta + \\alpha_{Richness} + \\beta^{Bias}$",
           "type*Richness" = "$\\beta + \\alpha_{Richness} + \\beta^{Bias} + \\alpha^{Bias}_{Richness}$",
           "type + ShannonEvenness" = "$\\beta + \\alpha_{Shannon\\,Evenness} + \\beta^{Bias}$",
           "type*ShannonEvenness" = "$\\left(\\overset{\\;\\,\\beta\\,\\;\\;+ \\alpha_{Shannon\\,Evenness}}{+ \\beta^{Bias} + \\alpha^{Bias}_{Shannon\\,Evenness}}\\right)$",
           "type*(Richness + ShannonEvenness)" = "$\\left(\\overset{\\;\\,\\beta\\,\\;\\;+ \\alpha_{Richness} + \\alpha_{Shannon\\,Evenness}}{+ \\beta^{Bias} + \\alpha^{Bias}_{Richness} + \\alpha^{Bias}_{Shannon\\,Evenness}}\\right)$") %>% latex2exp::TeX(out = "character")
  })) %>%
  mutate(form = factor(form,levels=unique(form)[order(nchar(unique(form)))])) %>% 
  rename(p.value = "Pr(>Chisq)") 

library(scales)
negative_log_trans <- function(base = exp(1)) {
  trans <- function(x) -log(abs(x), base)
  inv <- function(x) -base^(abs(x))
  trans_new(paste0("negative_log-", format(base)), trans, inv, 
            function(x,n) -log_breaks(n,base = base)(abs(x)), 
            domain = c(-Inf,-1e-100))
}

finalModelAIC <- finalModelAIC_dat %>% 
  mutate(Index = factor(Index)) %>% 
  group_by(Index) %>% 
  mutate(best = AIC == min(AIC)) %>% 
  mutate(AIC = ifelse(Index == "FEve" & model == 1,NA,AIC)) %>% 
  mutate(ann = ifelse(row_number()==1,paste0("(",letters[Index],")"),NA),
         ann_x = max(AIC,na.rm=T),
         form = as.numeric(form)) %>% 
  ggplot(aes(AIC,form,color=Index,group=Index)) +
  geom_path(size = 1) +
  geom_point(aes(fill=best,shape=best),
             size = 5, stroke = 1.5, shape = 21) +
  geom_text(inherit.aes = F,
            aes(ann_x,max(form),label = ann),
            position = position_nudge(y = -.25),
            family = "CMU Serif",
            fontface = "bold",
            size = 6) +
  scale_fill_manual(values = c("white", "black")) +
  # scale_shape_manual(values = c(16,21)) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(expand = expansion(0,.25),
                     breaks = 1:length(levels(finalModelAIC_dat$form)),
                     labels = parse(text = levels(finalModelAIC_dat$form)),
                     position = "right",
                     sec.axis = dup_axis()) + 
  scale_x_continuous(n.breaks = 4,
                     expand = expansion(.15,0)) +
  facet_wrap(~Index,scales = "free_x",nrow=1) +
  labs(y = "<span> Model Complexity <p style='font-size: 24px;'> &rarr; </p></span>", color=NULL, x = "Akaike's Information Criterion") +
  ggpubr::theme_pubr() +
  theme(text = element_text(family = "CMU Serif"),
        title = element_text(family = "CMU Serif",
                             face = "bold"),
        strip.text = element_text(family = "CMU Serif",
                                  face = "bold",
                                  size = 12),
        plot.margin = margin(.25,.25,.25,.25,"lines"),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   family = "CMU Serif"),
        axis.text.y = element_text(hjust = 0),
        axis.title.y.left = element_markdown(),
        axis.title.y.right = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        plot.caption = element_text(face = "italic",
                                    size = 9,
                                    hjust = 1),
        panel.spacing.x = unit(1.25,"lines"),
        legend.position = "none",
        aspect.ratio = 1.35)

ggsave("finalModelAIC.png",finalModelAIC,
       type = "windows", dpi = 500, width = 5, height = 2.05, scale = 2.25)

