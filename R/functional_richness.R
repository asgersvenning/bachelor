#' Calculates functional richness
#' 
#' Calculates functional richness from a species-trait matrix and possibly an abundance vector.
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param a optional numeric vector. Species-abundances.
#' @param scaleRich a logical. Are functional axes to be scaled?
#' @param absolute a logical. Should functional richness be standardized by maximum volume?
#' @return a number.
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#' 
#' @section Details:
#' When using 'absolute=T', the maximum volume is computed as the product of the functional trait axis ranges, and standardization is carried out by dividing the convex hull volume by this hypercube volume.
#'
#' @importFrom Rdpack reprompt
#' @importFrom geometry convhulln
#' @export

functional_richness <- function(x, a = rep(1, nrow(x)), scaleRich = F, absolute = F) {
  x <- x[a>0,] # Subset present species
  if (!absolute) x <- apply(x, 2, scale) # Scale and center values, such that richness is invariant to trait scales
  
  out <- geometry::convhulln(x,output.options = "FA")$vol / # Volume of convex hull
    ifelse(scaleRich, prod(abs(diff(apply(x,2,range)))),1)  # Possibly standardized by the volume of the minimum encapsulating (hyper-)cube
  
  return(out)
}
