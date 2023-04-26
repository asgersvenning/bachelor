#' Calculates functional divergence
#' 
#' Calculates functional divergence from a species-trait matrix and possibly an abundance vector. Weights should only be supplied, when using 'gower'=TRUE.
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param w numeric vector. A vector of length equal to columns in 'x', which specifies the variable weights. If missing, weights are equal.
#' @param a optional numeric vector. Species-abundances.
#' @param ch optional numeric or logical. Should the center of be based on the convex hull of the species, if true the convex hull is estimated on a 2-D PCoA. If a numeric the convex hull is estimateds on a N-D PCoA.
#' @param gower a logical. Calculate entropy based on Gower dissimilarity as opposed to euclidean distance.
#' @param detailed a logical. If true the partial divergences and convex hull vertices are returned.
#' @return a number. Or if 'detailed'=TRUE, a list containing; the community functional divergences, the partial divergences and a boolean indicating whether a species is in the set of convex hull vertices for all present species in each site (missing values are missing species).
#' 
#' @section Details:
#' This functions implements functional divergence as it is defined in \insertCite{Villeger2008}{asgerbachelor}.
#' It handles both euclidean distance and Gower dissimilarity (however not with ordered factors or asymmetric binary variables).
#' 
#' The argument 'ch' defines a maximum trait dimensionality, which will be lowered on a site by site basis, if the number of species is lower than the number of (ordinated) traits.
#' OBS: If 'ch'=FALSE, then no convex hull will be calculated, and instead all site-present species will be used for center of mass calculation.
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#' 
#' @importFrom magrittr %>% 
#' @importFrom Rdpack reprompt
#' @importFrom grDevices chull
#' @importFrom stats cmdscale
#' @importFrom stats dist
#' @export

functional_divergence <- function(x, w, a = rep(1, nrow(x)), ch = T, gower = T, detailed = F) {
  # Attempt to coerce x and a to matrices.
  if (is.vector(a)) a <- matrix(a, nrow = 1)
  if (inherits(x,"data.frame")) x <- as.matrix(x)
  if (inherits(a,"data.frame")) a <- as.matrix(a)
  
  if (!is.matrix(a)) stop("Unable to coerce 'a' to matrix.")
  if (!is.matrix(x)) stop("Unable to coerce 'x' to matrix.")
  
  if (missing(w)) w <- rep(1,ncol(x))
  if (!is.vector(w) | !is.numeric(w)) stop("'w' must be a numeric vector.")
  
  if (length(ch) != 1) stop("'ch' must have length 1.")
  ch <- if (is.numeric(ch)) ch else if (isTRUE(ch)) 2 else if (isFALSE(ch)) 0 else stop("'ch' must be a logical or numeric of length 1.")
  
  # Missing abundances are replaced by 0, 
  a <- replace(a, is.na(a), 0)
  # and all abundances are turned into relative abundances.
  a <- a / rowSums(a)
  
  # Pairwise dissimilarities between all species. Will be used later.
  d <- if (gower) gowerDissimilarity(x, w) else as.matrix(dist(x))
  
  # Calculates a "convex hull" matrix, where each rows are communities and columns are species, 
  # while entries are the boolean values for if a species is a vertex in the convex hull. 
  chMat <- t(apply(a, 1, function(com) {
    # Find the subset of species in the current community
    comSP <- which(com != 0)
    # Extract their precomputed pairwise dissimilarity ...
    comD <- d[comSP,comSP]
    # and trait matrix.
    comX <- x[comSP,]
    
    # Warn the user if a community contains a single species
    if (length(comSP) %in% 1:2) {
      warning("Community with 1 or 2 single species.")
      return(if (length(comSP) == 1) comX else colMeans(comX))
    } 
    
    # Calculate the indices for the convex hull over the set of species in the current community.
    chInd <- if (ch != 0) { 
      # Choose the appropriate coordinates to calculate the convex hull on, 
      # given the number of species, traits and the value of ch:
      # First if the number of traits is lower than or equal to the desired dimensionality of the convex hull ... 
      coord <- if (ncol(x) <= ch) {
        # and the number of species (points) is greater than the number of traits (dimensions),
        # no dimensionality reduction is performed.
        if (nrow(comX) > ncol(comX)) {
          comX
        } else {
          # If the number of traits is greater than the number of species, 
          # then the traits are ordinated to 1 less dimension than there are species.
          cmdscale(comD, nrow(comX)-1)
        }
      } else {
        # And if there are more traits, than the desired dimensionality of the convex hull,
        # then the traits are ordinated to the desired dimensionality.
        cmdscale(comD, ch)
      }
      
      # Calculate the convex hull using the method appropriate for the number of dimensions.
      if (ncol(coord) == 1) {
        # For 1 dimension, the "convex hull" simply consists of the maximum and minimum points.
        c(which.min(coord), which.max(coord)) 
      } else if (ncol(coord) == 2) {
        # For 2 dimensions, the base R function chull can be used, which is fast and effecient.
        chull(coord)
      } else {
        # For N dimensions, the convhulln function from the geometry package must be used.
        unique(as.vector(geometry::convhulln(coord)))
      }
    } else {
      # The user is also provided the option of not calculating the convex hull, but instead simply using all species.
      1:nrow(comX)
    }
    
    out <- rep(0,length(com))
    
    out[comSP][chInd] <- 1
    
    out
  }))
  
  centers <- (chMat / rowSums(chMat)) %*% x 
  
  # Calculate all distances to center of gravity
  dcg <- if (gower) {
    gowerDissimilarity(x, w, centers)
  }
  else {
    t(apply(centers, 1, function(z) sqrt(colSums((t(x) - z)^2))))
  }
  
  # Mean distance to center of gravity.
  mdcg <- rowMeans(dcg * ifelse(a==0,NA,1), na.rm=T)
  
  # Distance anomaly from center of gravity.
  ddcg <- dcg - mdcg
  
  # Weighted sum of distance anomalies.
  wddcg <- rowSums(a * ddcg) 
  
  # Weighted sum of absolute distance anomalies.
  awddcg <- rowSums(a * abs(ddcg))  
  
  # Finally the divergence is calculated using its defining formula.
  div <- (wddcg + mdcg)/(awddcg + mdcg)
  
  # The output type is defined by the user argument 'detailed'.
  out <- if (!detailed) {
    div
  } else {
    chMat[which(a==0)] <- NA
    list(FDiv = div, partial = ddcg, convexhulls = chMat)
  }
  
  return(out)
}
