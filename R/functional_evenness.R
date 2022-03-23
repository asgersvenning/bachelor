#' Calculates functional evenness
#' 
#' Calculates functional evenness from a species-trait matrix and possibly an abundance vector.
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param w numeric vector. A relative variable importance vector. Not used in case 'gower' = FALSE.
#' @param a optional numeric vector. Species-abundances.
#' @param gower a logical. Calculate entropy based on Gower dissimilarity as opposed to euclidean distance.
#' @return a number. A double between 0 and 1.
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#' 
#' @section Details:
#' The minimum spanning tree is constructed from a pairwise similiarity matrix, which is an inefficient way to do so.
#' It seems however, that efficient implementations of minimum spanning trees escape the algorithmic landscape of R.
#' 
#' Species with an abundance of 0, are automatically filtered, before any computations.
#' 
#' @importFrom magrittr %>% 
#' @importFrom Rdpack reprompt
#' @importFrom spdep knearneigh
#' @importFrom spdep knn2nb
#' @importFrom spdep nb2listw
#' @importFrom ape mst
#' @export

functional_evenness <- function(x, w, a = rep(1, nrow(x)), gower = T) {
  x <- x[a>0,] # Subset present species
  a <- a[a>0]  # Subset present species
  a <- a/sum(a)# Relative abundances
  
  S <- length(a)# Number of species
  Ed <- 1/(S-1) # Expected edge length in minimum spanning tree
  
  dx <- if (gower) FD::gowdis(x,w) else dist(x)
  l <- unclass(ape::mst(dx))==1 # Logical encoding of minimum spanning tree
  li <- which(l, arr.ind = T) # Indices of mst-vertices
  li.lower <- li[diff(t(li))<0,] # Subset lower triangle
  
  EW <- apply(li.lower, 1, function(y) dx[y[1],y[2]]/sum(a[y]))
  
  PEW <- EW/sum(EW) # Relative weighted length of mst-edge
  
  out <- (sum(pmin(PEW,Ed))-Ed)/(1-Ed) 
  
  return(out)
}