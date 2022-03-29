#' Calculates functional dispersion
#' 
#' Calculates functional dispersion from a species-trait matrix and possibly an abundance vector
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param w numeric vector. A vector of length equal to columns in 'x', which specifies the variable weights. If missing, weights are equal.
#' @param a optional numeric vector. Species-abundances.
#' @param gower a logical. Calculate entropy based on Gower dissimilarity as opposed to euclidean distance.
#' @param returnPartial an optional logical. If true the distances of each species to each center is returned along with the functional dispersion. 
#' @return a number.
#' 
#' @section Details:
#' This functions implements functional dispersion, as the mean distance/dissimilarity to the community weighted mean (traits). 
#' 
#' It also supplies the option of returning the "partial" dispersions, i.e. the distance/dissimilarity of each species in each site, from that center.
#' 
#' It also supplies the option of calculating the center of mass based on the convex hull of the first 2 or N principal coordinates. Caution should be used, when using this option, since this deviates from the definition of functional dispersion.
#' 
#' @importFrom magrittr %>% 
#' @importFrom Rdpack reprompt
#' @export

functional_dispersion <- function(x, w, a = rep(1, nrow(x)), ch = F, gower = T, returnPartial = F) {
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
  center <- if (!isFALSE(ch)) {
    if (isTRUE(ch) | ch==2) {
      x[chull(cmdscale({if (gower) FD::gowdis(x,w) else dist(x)})),]
    } else if (is.numeric(ch)) {
      x[geometry::convhulln(cmdscale({if (gower) FD::gowdis(x,w) else dist(x)},ch), output.options = "FA")$hull %>% 
          as.vector %>% 
          unique,]
    } else {
      stop("'ch' must be a logical or numeric.") 
    }
  } else {
    x
  }
  
  # Calculates the mean
  center <- as.matrix(center) %>% 
    inset(is.na(.),0) %>% 
    t %>% 
    multiply_by_matrix(t(a))  %>% 
    t
  
  # Calculate all distances to center of gravity
  dcg <- if (gower) gowerDissimilarity(x,w,colMeans(center,na.rm=T)) else sqrt(colSums((t(x) - colMeans(center,na.rm=T))^2))
  
  # Mean distance to center of gravity
  mdcg <- as.vector((dcg %*% t((a>0) %>% {./rowSums(.)}))) 
  
  out <- if (!returnPartial) mdcg else list(FDiv = mdcg, partial = dcg)
  
  return(out)
}
