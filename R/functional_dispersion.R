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

functional_dispersion <- function (x, w, a = rep(1, nrow(x)), ch = F, gower = T, returnPartial = F) 
{
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
  
  d <- if (gower) gowerDissimilarity(x, w) else dist(x)
  
  centerInd <- if (!isFALSE(ch)) {
    if (isTRUE(ch) | ch == 2) {
      chull(cmdscale(d))
    }
    else if (is.numeric(ch)) {
      unique(as.vector(geometry::convhulln(cmdscale(d, ch))))
    }
    else {
      stop("'ch' must be a logical or numeric.")
    }
  }
  else {
    1:nrow(x)
  }
  center <- a[,centerInd] %*% as.matrix(x[centerInd,])
  
  dcg <- if (gower) {
    gowerDissimilarity(x, w, center)
  }
  else {
    t(apply(center, 1, function(z) sqrt(colSums((t(x) - z)^2))))
  }
  
  partialDisp <- ifelse(a==0,NA,a) * dcg
  
  dispersion <- rowSums(partialDisp, na.rm=T)
  
  out <- if (!returnPartial) 
    dispersion
  else {
    partialDisp[partialDisp==0] <- NA
    list(FDis = dispersion, partial = partialDisp)
  }
  return(out)
}
