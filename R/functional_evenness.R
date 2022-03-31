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

functional_evenness <- function (x, w, a = rep(1, nrow(x)), gower = T) {
  if (is.vector(a)) 
    a <- matrix(a, nrow = 1)
  if (inherits(x, "data.frame")) 
    x <- as.matrix(x)
  if (inherits(a, "data.frame")) 
    a <- as.matrix(a)
  if (!is.matrix(a)) 
    stop("Unable to coerce 'a' to matrix.")
  if (!is.matrix(x)) 
    stop("Unable to coerce 'x' to matrix.")
  if (missing(w)) 
    w <- rep(1, ncol(x))
  if (!is.vector(w) | !is.numeric(w)) 
    stop("'w' must be a numeric vector.")
  a <- replace(a, is.na(a), 0)
  a <- a/rowSums(a)
  
  dx <- as.matrix({if (gower) gowerDissimilarity(x, w) else dist(x)})
  
  apply(a, 1, function(com) {
    S_list <- which(com>0)
    S <- length(S_list)
    Ed <- 1/(S - 1)
    
    aCom <- com[S_list]
    dCom <-dx[S_list,S_list]
    
    l <- unclass(ape::mst(dCom)) == 1
    li <- which(l, arr.ind = T)
    li.lower <- li[diff(t(li)) < 0, ]
    EW <- apply(li.lower, 1, function(y) dCom[y[1], y[2]]/sum(aCom[y]))
    PEW <- EW/sum(EW)
    
    
    (sum(pmin(PEW, Ed)) - Ed)/(1 - Ed)
  })
}