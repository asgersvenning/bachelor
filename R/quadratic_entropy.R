#' Calculates functional identity
#' 
#' Calculates functional identity from a species-trait matrix and possibly an abundance vector.
#' 
#' @param x numeric matrix. Species-trait matrix.
#' @param a optional numeric vector. Species-abundances.
#' @param raoScale a logical. Is quadratic entropy divided by maximum value (ade4::divcmax)?
#' @param Botta_Dukat a logical. Are distances standardized by number of trait-axis (dimensions)?
#' @param approxRao an integer. If NOT 0, number of species pairs to estimate quadratic entropy on.
#' @return a number.
#' 
#' @references
#' \insertRef{Villeger2008}{asgerbachelor}
#'  
#'  
#' @importFrom magrittr %>% 
#' @importFrom Rdpack reprompt
#' @export

quadratic_entropy <- function(x, a = rep(1, nrow(x)), raoScale=T, Botta_Dukat=F, approxRao=0) {
  if (raoScale & Botta_Dukat) stop("Only one type of scaling should be applied.") # Either standardize by number of traits, or maximum raoQ, not both.
  if (approxRao!=0) {
    if (raoScale) stop("Unable to use divcmax on less than full dissimilarity matrix")
    return(sapply(1:approxRao, function(y) { # Sample (n=approxRao) pairs weighted by abundance, and calculate half-mean
      ind <- sample.int(nrow(x),2,T, prob = n) # Pair indices
      sum((x[ind[1],]-x[ind[2],])^2) # Pair dissimilarity
    }) %>% 
      mean / ifelse(Botta_Dukat,ncol(x),1) / 2) # Mean pair dissimilarity divided by 2 (and possible number of traits)
  }
  
  dx <- dist(x)  # Dissimilarity matrix
  
  D <- as.matrix(dx)^2/ifelse(Botta_Dukat,ncol(x),1) # Possibly standardized by number of traits
  
  out <- sum(D * outer(n,n)) / 2 / sum(n)^2 / # Abundance weighted sum of lower triangle of dissimilarity matrix
    ifelse(raoScale, ade4::divcmax(dx)$value,1) # Possible scaled by maximum possible raoQ
  
  return(out)
}