train <- ecoCommunitiesOrd %>% 
  group_by(eco) %>% 
  transmute(
    inTrain = sample(c(T,F),n(),T,prob = c(.8,.2))
  ) %>% 
  pull(inTrain)

SVM <- ksvm(eco ~ ., data = ecoCommunitiesOrd[train,],
            kernel = "laplacedot",
            scaled = F,
            C = 30)

SVM_novel <- ksvm(eco ~ ., data = ecoCommunitiesOrd,
            kernel = "laplacedot",
            scaled = F,
            nu = .1,
            type = "one-svc")

tEval <- ecoCommunitiesOrd[,2:19] %>% 
  summarize(across(everything(), ~rnorm(1000,mean(.x),sd(.x)))) %>% 
  mutate(isReal = predict(SVM_novel,.) %>% as.vector) %>% 
  filter(isReal) %>% 
  select(!isReal) %>% 
  mutate(class = predict(SVM,.))


tEval %>%
  group_by(int1 = cut_interval(V1,16,boundary = 0, dig.lab = 10),
           int2 = cut_interval(V2,16,boundary = 0, dig.lab = 10)) %>% 
  summarize(
    maxclass = names(which.max(table(class))),
    prob = mean(class==maxclass),
    .groups = "drop"
  ) %>% 
  mutate(across(c(int1,int2), ~str_extract_all(.x,"[[:digit:]\\.\\-]+") %>% 
                  unlist %>% 
                  as.numeric %>% 
                  matrix(nrow = 2) %>% 
                  colMeans)) %>% 
  ggplot(aes(int1,int2,fill=maxclass)) +
  geom_raster() +
  scale_alpha_continuous(range = 0:1,
                         trans = "logit")


tEval %>% 
  ggplot(aes(V1,V4,group = -1L)) +
  geom_point() + 
  ggforce::geom_voronoi_tile(aes(fill=class)) +
  coord_cartesian(expand = F, clip = "off")
