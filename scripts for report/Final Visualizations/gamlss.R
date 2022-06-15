library(tidyverse)
library(extrafont)
library(nlme)
library(multcomp)
library(rstatix)
library(gamlss)
library(patchwork)
library(emmeans)
library(mgcv)
library(magrittr)

dat <- tibble(mu = sample(0:10 * 5,10,T),
              sd = runif(10,.5,5))

nObs <- 50

sampleDat <- dat %>% 
  group_by(grp = row_number() %>% factor,
           pred = sample(c(-20:10),n(),T)) %>% 
  summarize(pred = runif(nObs,pred/2,5+pred/2),
            sample = 0.25*pred^2 - pred - 5 * sin(pred + 5) + rnorm(nObs,mu,sd),
            ground = mu) %>% 
  ungroup %>% 
  mutate(grp = factor(grp))

list(ggplot(sampleDat, aes(pred,sample,color=grp)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x)),
  ggplot(sampleDat,aes(grp,sample)) +
    geom_violin() +
    stat_summary(fun.data = mean_sdl)) %>% 
  wrap_plots(ncol = 1)


# if (exists("tGAM")) rm(tGAM)
# for (i in 1:20) {
#   res <- if (exists("tGAM")) tibble(res = residuals(tGAM),
#                                     grp = sampleDat$grp) %>%
#     group_by(grp) %>%
#     mutate(var = var(res)) %>%
#     pull(var) else rep(1, nrow(sampleDat))
# 
#   tGAM <- gam(sample ~ s(pred) + grp,
#                 data = sampleDat,
#                 weights = 1 / res
#   )
# 
#   if (i == 5) {
#     ftGAM <- tGAM
#     rm(tGAM)
#   }
# }

ftGAM <- gam(list(sample ~ s(pred) + grp,
                  ~ grp),
            data = sampleDat,
            family = gaulss())

# tibble(res = residuals(ftGAM,"pearson"),
#        fit = fitted(ftGAM),
#        grp = sampleDat$grp) %>%
#   ggplot(aes(grp,res)) +
#   geom_violin() +
#   stat_summary(fun.data = mean_sdl,
#                color = "firebrick") +
#   baseTheme +
#   theme(aspect.ratio = 1)

nGroup <- sum(str_detect(names(coef(ftGAM)),"grp[:digit:]+$")) + 1
nVar <- length(coef(ftGAM))
tMAT <- matrix(0,(sum(1:nGroup) - nGroup),nVar) %>% 
  set_colnames(names(coef(ftGAM)))
tMAT[matrix(c(1:nrow(tMAT),rep(1:nGroup,(nGroup-1):0)),ncol=2)] <- -1
tMAT[matrix(c(1:nrow(tMAT),unlist(sapply(1:(nGroup-1), function(x) (x+1):nGroup))),ncol=2)] <- 1
rownames(tMAT) <- tMAT %>% apply(1,function(x) paste0(names(which(x != 0)), collapse = " - "))
tMAT[,1] <- 0


summary(glht(ftGAM,tMAT),test = adjusted("BH"))$test[3:6] %>% 
  as_tibble %>% 
  mutate(comp = names(coefficients)) %>% 
  filter(pvalues > .05) 
  
sampleDat %>% 
  group_by(ground) %>% 
  summarize(
    groups = paste0(unique(grp) %>% as.numeric %>% sort,collapse=", ")
  ) %>% 
  arrange(as.numeric(str_extract(groups,"[:digit:]+"))) %>% 
  filter(str_detect(groups,","))

