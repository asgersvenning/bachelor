#' Calculates functional divergence
#' 
#' Calculates functional divergence from a species-trait matrix and possibly an abundance vector
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param a optional numeric vector. Species-abundances.
#' @return a number.
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#' 
#' @importFrom magrittr %>% 
#' @importFrom Rdpack reprompt
#' @export

functional_divergence <- function(x, a = rep(1, nrow(x))) {
  x <- x[n>0,] # Subset present species
  n <- n[n>0]  # Subset present species
  n <- n/sum(n)# Calculate relative abundances
  
  l <- chull(x) # Calculate indices of convex hull
  
  # cg <- colMeans(x[l,]) # Calculate center of gravity using the convex hull
  dcg <- apply(x,2,function(y) (y-mean(y[l])^2)) %>% 
    rowSums %>% 
    sqrt # Calculate all distances to center of gravity
  mdcg <- mean(dcg) # Mean distance to center of gravity
  ddcg <- dcg - mdcg # Distance anomaly from center of gravity
  
  wddcg <- sum(ddcg*n) # Weighted sum of distance anomalies
  awddcg <- sum(abs(ddcg)*n)  # Weighted sum of absolute distance anomalies
  out <- (wddcg + mdcg)/(awddcg + mdcg)
  
  return(out)
}