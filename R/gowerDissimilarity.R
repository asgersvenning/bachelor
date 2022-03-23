#' Calculates Gower Dissimilarity either between observations or sets.
#' 
#' This function is slightly slower than FD::gowdis, however instead it provides the flexibility of calculating Gower Dissimilarity between sets of observations.
#' An obvious application of this function, is computing the distribution of distances to the center of mass. 
#' 
#' @param set numeric matrix. Species-trait matrix.
#' @param w numeric vector. Relative trait weights, can be missing in which case all weights are equal.
#' @param set2 numeric matrix or vector. Species-trait matrix of the comparison set, can be missing in which case it is equal to 'set'.
#' @return A matrix of the pairwise Gower Dissimilarities between species (rows) in the first and second set. 
#' If the second set is missing, simply a pairwise Gower Dissimilarity matrix.
#' If the second set contains only a single "species", a vector of Gower Dissimilarities between species in the first set and the single "species" in the second set.
#' 
#' @section Details:
#' This functions removes missing values from the Gower Dissimilarity calculation and reweights the trait, based on the sum of the weights of non-missing traits.
#' This step is important to ensure that the calculation can be carried out, even if there are missing values scattered in the trait matrix, however it can also be slightly misleading, since the calculation can be carried out even if all but a single trait is missing, in which case the dissimilarity is based solely on the remaining trait. 
#' 
#' @importFrom magrittr %>% 
#' @importFrom magrittr set_colnames
#' @importFrom magrittr set_names
#' @importFrom magrittr inset
#' @importFrom magrittr multiply_by_matrix
#' @importFrom magrittr divide_by
#' @importFrom magrittr not
#' @importFrom tibble tibble
#' @importFrom tibble as_tibble
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr select
#' @importFrom purrr map
#' @importFrom tidyr unnest
#' @importFrom tidyr replace_na
#' @importFrom Rdpack reprompt
#' @export
#' 
#' @examples 
#' # Obviously "tidyverse" needs to be loaded for most of these examples to work.
#' createSyntheticMat <- function(x,n,na.prob=0.05) {
#'   vec <- runif(x*n)*sample(c(1,NA),x*n,T,prob=c(1-na.prob,na.prob))
#'   
#'   matrix(vec,ncol=x,nrow=n)
#' }

#' compGowerDis <- lapply(1:15, function(x) {
#'   lapply((1:25)^2, function(y) {
#'     dat <- createSyntheticMat(x,y) %>% 
#'       gower_traits()
#'     
#'     tibble(
#'       res = list(system.time(gowdis_spec(dat$traits,dat$weights)),
#'                  system.time(FD::gowdis(dat$traits,dat$weights))) %>% 
#'         lapply(as.list) %>% 
#'         lapply(as_tibble),
#'       type = c("Mine","FD"),
#'       variables = x,
#'       species = y)
#'   })
#' })
#' 
#' tibble(compGowerDis) %>% 
#'   unnest(1) %>% 
#'   unnest(1) %>% 
#'   unnest(res) %>%
#'   select(-c(4:5)) %>% 
#'   pivot_longer(1:3,names_to="timerType",values_to = "time") %>% 
#'   ggplot(aes(variables,species,z=time)) +
#'   stat_summary_2d(bins = 10) +
#'   scale_fill_viridis_c(option = "A",
#'                        trans = "log10") +
#'   facet_grid(rows = vars(type), cols = vars(timerType))
#' 
#' PLANTS_tG <- PLANTS_traits %>%
#' filter(rowMeans(across(!plants_code,is.na))<.5) %>%
#' gower_traits(T)
#' 
#' allGowdisMine <- gowerDissimilarity(as.matrix(PLANTS_tG$traits[,-1]),PLANTS_tG$weights)
#' 
#' minePCoA <- ape::pcoa(allGowdisMine,rn=PLANTS_tG$traits[[1]])
#' 
#' minePCoA$values[,3:4] %>% 
#'   as_tibble %>% 
#'   mutate(nPC = row_number()) %>% 
#'   ggplot(aes(nPC)) +
#'   geom_col(aes(y=Rel_corr_eig,fill=Rel_corr_eig>Broken_stick),
#'            width = 1,
#'            position = "dodge") +
#'   geom_path(aes(y=Broken_stick),color="red")
#' 
#' minePCoA %>% biplot()
#' 
#' distToCOM <- gowerDissimilarity(PLANTS_tG$traits[,-1],PLANTS_tG$weights,PLANTS_tG$traits[,-1] %>% colMeans(na.rm=T))
#' distToCOH <- gowerDissimilarity(PLANTS_tG$traits[,-1],PLANTS_tG$weights,PLANTS_tG$traits[chull(minePCoA$vectors[,1:2]),-1] %>% colMeans(na.rm=T))
#' 
#' par(mfrow = c(1,3))
#' distToCOM %>%
#'   hist(seq(0,.5,.01),main="Trait-Vector Gower Dissimilarity\nwith Center of Mass",xlim=c(0,.5),)
#' 
#' distToCOH %>%
#'   hist(seq(0,.5,.01),main="Trait-Vector Gower Dissimilarity\nwith Center of Convex Hull",xlim=c(0,.5),)
#' 
#' minePCoA$vectors %>%
#'   raise_to_power(2) %>%
#'   rowSums %>%
#'   sqrt %>%
#'   hist(seq(0,.5,.01),main = "PCoA-Vector Magnitudes",xlim=c(0,.5))


gowerDissimilarity <- function(set, w=rep(1,ncol(set)), set2) {
  # If set2 is missing, the function assumes that a pairwise comparison is desired.
  if (missing(set2)) set2 <- set
  
  # Attempt to coerce the sets to matrices
  # If either set or set2 is a vector, convert to a 1-row matrix
  if (is.vector(set)) set <- matrix(set,nrow=1)
  if (is.vector(set2)) set2 <- matrix(set2,nrow=1)

  if (inherits(set,"data.frame")) set <- as.matrix(set)
  if (inherits(set2,"data.frame")) set2 <- as.matrix(set2)
  
  if (!is.matrix(set)) stop("Unable to coerce 'set' to matrix.")
  if (!is.matrix(set2)) stop("Unable to coerce 'set2' to matrix.")
  
  if (!is.vector(w) | !is.numeric(w)) stop("'w' must be a numeric vector.")
  
  # Calculate relative weights
  w <- w %>% 
    # Missing weights are interpreted as 0.
    inset(is.na(.),0) %>% 
    divide_by(sum(.))
  
  # To reduce memory usage and keep the calculations to matrix multiplication, 
  # in fact the entire calculation could be carried out with a single tensor contraction operation.
  # However the size of the intermediate tensor would be equal to N(set1) x N(set2) x N(w).
  # This way calculation is converted to N(set) matrix multiplication, 
  # with intermediate matrices of size of N(set2) x N(w).
  apply(set, 1, function(x) {
    delta <- abs(t(set2) - x) %>% 
      t 
    
    delta %>% 
      inset(is.na(.),0) %>% 
      multiply_by_matrix(w) %>% 
      divide_by(delta %>% 
                  is.na %>% 
                  not %>% 
                  multiply_by_matrix(w))
  })
}


 
