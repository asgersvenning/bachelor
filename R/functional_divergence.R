#' Calculates functional divergence
#' 
#' Calculates functional divergence from a species-trait matrix and possibly an abundance vector
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param w numeric vector. A vector of length equal to columns in 'x', which specifies the variable weights. If missing, weights are equal.
#' @param a optional numeric vector. Species-abundances.
#' @param gower a logical. Calculate entropy based on Gower dissimilarity as opposed to euclidean distance.
#' @return a number.
#' 
#' @section Details:
#' This functions implements functional divergence as it is defined in \insertCite{Villeger2008}{asgerbachelor}, 
#' based on a weighted Gower dissimilarity measure (which is equivalent to the calculation in FD::gowdis, expect it doesn't handle ordered factors or assymetric binary variables).
#' 
#' It also supplies the option of calculating the center of mass based on the convex hull of the first two principal coordinates.
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#' 
#' @importFrom magrittr %>% 
#' @importFrom Rdpack reprompt
#' @export

functional_divergence <- function(x, w, a = rep(1, nrow(x)), ch = F, gower = T) {
  # Attempt to coerce x and a to matrices.
  if (is.vector(a)) a <- matrix(a, nrow = 1)
  if (inherits(x,"data.frame")) x <- as.matrix(x)
  if (inherits(a,"data.frame")) a <- as.matrix(a)
  
  if (!is.matrix(a)) stop("Unable to coerce 'a' to matrix.")
  if (!is.matrix(x)) stop("Unable to coerce 'x' to matrix.")
  
  if (missing(w)) w <- rep(1,ncol(x))
  if (!is.vector(w) | !is.numeric(w)) stop("'w' must be a numeric vector.")
  
  a <- replace(a, is.na(a), 0)
  a <- a / rowSums(a)
  
  # Calculates which point to use for center of mass calculation. 
  center <- if (ch) x[chull(cdmscale(FD::gowdis(x,w))),] else x
  
  # Calculate all distances to center of gravity
  dcg <- if (gower) gowerDissimilarity(x,w,colMeans(center,na.rm=T)) else sqrt(colSums((t(x) - colMeans(center,na.rm=T))^2))
  
  mdcg <- rowMeans(dcg * (a>0)) # Mean distance to center of gravity
  ddcg <- dcg - mdcg # Distance anomaly from center of gravity
  
  wddcg <- rowSums(a %*% ddcg) # Weighted sum of distance anomalies
  awddcg <- rowSums(abs(ddcg)*a)  # Weighted sum of absolute distance anomalies
  out <- (wddcg + mdcg)/(awddcg + mdcg)
  
  return(out)
}
